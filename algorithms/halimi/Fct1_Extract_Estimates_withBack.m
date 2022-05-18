function [Dscales, Rscales, ProfileT , ProfileN, sumB,Depth_CM] = Fct1_Extract_Estimates_withBack(Y,I_resol,F,IRFw,Tbin,Attack  ,trailing,mask)
%% input: I_resol, Y, F, IRFw,
%% output: Dscales, Rscales, Back

[row,col,T,L] = size(Y);
N = row*col;
R = length(I_resol);
Yvectorized    = permute(reshape(Y, [row*col,T,L]), [2,1,3]); % T N L
if(T > 800 )
    IRFw = 0;
else
    Yvectorized   = cat(1,Yvectorized(IRFw:-1:1,:,:),Yvectorized,Yvectorized(end-IRFw+1:end,:,:) );
end;

%%%% preparation of IRF
if(prod(size(F))> (T*L))
    h = squeeze(F(:,floor(T/2),:));
else
    h = F;
end

if size(h,1)>(T+2*IRFw)
    h = h(1:sY3,:);
else
    d = zeros((T+2*IRFw),size(h,2));
    d(1:size(h,1),:) = h;
    h = d;
end
[~,attack] = max(h);
for ell=1:L
    h(:,ell) = circshift(h(:,ell),[-attack(ell),0]);
end
h      = flipud(h);
Hf     = fft(h); % T  1

%% Convolution in time
% % tic
% % Yf     = fft(Yvectorized);%; zeros(sum(h>0.01),size(Ymat,2)) ]  ); % T  N
% % for ell=1:L
% %     YHf(:,:,ell) =  Yf(:,:,ell).*Hf(:,ell);
% % end
% % Z_F         = ifft(YHf,'symmetric'); % T  N
% % toc
% % % keyboard
% % tic
Z_F = zeros(T+2*IRFw,N,L);
for ell=1:L
    Hf     = fft(h(:,ell)); % T  1
    Z_F(:,:,ell) = ifft( fft(Yvectorized(:,:,ell)) .*Hf  ,'symmetric');
end

clear Yf Yvectorized Hf Y
Z_F         = permute(Z_F,[2,1,3]); % N T L
Z_Fextended = reshape(Z_F,[row,col,size(Z_F,2),L]);
clear Z_F


%% Convolution in space 
for ell=1:L % wavelengths
    for r=1:R % multiscale
        %         r
        if(I_resol(r)>1)
            %             Z_Fextended(:,:,:,ell,r)  = imboxfilt3(Z_Fextended(:,:,:,ell,1) , [I_resol(r) I_resol(r) 1],'Padding' ,'replicate');
            Z_Fextended(:,:,:,ell,r)   = convn(Z_Fextended(:,:,:,ell,1) , ones(I_resol(r),I_resol(r))/I_resol(r)^2,'same');
            mask2 =  conv2(mask(:,:,ell,1) , ones(I_resol(r),I_resol(r))/I_resol(r)^2,'same');
            mask2(mask2==0) =1;
            Z_Fextended(:,:,:,ell,r)   =  Z_Fextended(:,:,:,ell,r)  ./mask2;
        end
    end
end
% max(abs(Z_Fextended2(:) - Z_Fextended3(:)))
% keyboard

if(IRFw>0)
    Z_Fextended       = Z_Fextended(:,:,IRFw+1:T+IRFw,:,:); % r c T  L R
end
Z_FextendedReshape = reshape(Z_Fextended,[row*col,T,L,R]); % N T  L R
clear Z_Fextended

%% Estimate background
% if(~exist('Back'))
Ysort     = sort(Z_FextendedReshape(:,:,:,end)); % N  T L
ProfileT  = median(Ysort(1:floor(0.2*N),:,:)); % 1 T L
ProfileT0 = ProfileT;
ProfileT  = (ProfileT-mean(ProfileT,2)); % 1 T L
ProfileN  = median(Z_FextendedReshape(:,1:T,:,end),2); %N 1 L
%     Back      = max(0,ProfileT + ProfileN); %N T L
% end

%% Estimate depths
sumB = sum(sum(sum( max(0,ProfileT + ProfileN))));
Z_Cleaned  = max(Z_FextendedReshape -  max(0,ProfileT + ProfileN), 0); % N T L R
clear Z_FextendedReshape
[v, Depth] = max(sum(Z_Cleaned,3),[],2); % N 1 1 R
Depth      = squeeze(Depth); % N R

%% Estimate Intensities
% Yvectorized=permute(Yvectorized,[2,1,3]);% N  T L
% rr = (sum(Yvectorized,2)); % rr
% Back2=reshape(Back,[row,col,T,L]);
%  for n = 1:N
% DerL = 1-Yvectorized(n,:,:)./ (rr(n,1,:) + Back(n,:,:)./F(:,Depth(n,r),:))
% rr(n,1,:) = rr(n,1,:) - DerL;
%  end

if(T<801)
    for r=1:R % multiscale
        %%% Compute intensity
        for n = 1:N
            %     int = max(1,Depth(n)+minh):min(T,Depth(n)+maxh);
            %         int        = max(1,Depth(n,r) - 7):min(T,Depth(n,r)+40);
            for ell=1:L % Wavelength
                int             = max(1,Depth(n,r) -trailing(ell)):min(T,Depth(n,r)+Attack(ell));
                weights         = Z_Cleaned(n,int,ell,r);%Z_F(int,n,:,r);%    max(0,Zdump0(int,n) - bb(:) -  sqrt(ddd)*1*sqrt(bb(:)));
                Intens(n,ell,r) = sum(weights,2);%sum(weights,1);
                Depth_CM(n,r)   = sum(int.*weights)/max(Intens(n,ell,r),eps) ;
            end
        end
    end
    Rscales = reshape(Intens,[row,col,L,R]);
else
    Rscales  = reshape(sum(Z_Cleaned,2),[row,col,L,R]); % r c  L R
    Depth_CM = Depth;
end

% keyboard
clear Z_Cleaned
Dscales  = reshape(Depth,[row,col,R])* Tbin*3*10^8/2;
Depth_CM  = reshape(Depth_CM,[row,col,R])* Tbin*3*10^8/2;

%
% figure;for i=1:3 subplot(2,3,i);imagesc(Dscales(:,:,i) );end
% figure;for i=1:3 subplot(2,3,i);imagesc(Rscales(:,:,i) );end
% Ym = reshape(Y(:,:,:,1),N,T);
% figure;plot(sum(Ym));hold on;plot(squeeze(sum(Back,1)))
% keyboard

