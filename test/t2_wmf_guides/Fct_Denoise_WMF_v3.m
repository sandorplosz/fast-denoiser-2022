function  [Dscales, Rscales, sumB, D_uncert, R_uncert, W, Dep, Int] = ...
    Fct_Denoise_WMF_v3(Y, Var2,I_resol,F,IRFw,Tbin,Attack  ,trailing,...
        Neighbours,indGraph,local_shifts,Dtype,ThreshDep,GuideD,GuideI,CvPC, usewmf, useguide)

switch Dtype % check dta format
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Histogram'
        %% input: I_resol, Y, F, IRFw,
        %% output: Dscales, Rscales, Back
        
        [row,col,T,L] = size(Y);
        N = row*col;
        R = length(I_resol);
        
        %% Expand Y
        Yvectorized    = permute(reshape(Y, [row*col,T,L]), [2,1,3]); % T N L
        if(T > 800 )
            IRFw = 0;
        else
            Yvectorized   = cat(1,Yvectorized(IRFw:-1:1,:,:),Yvectorized,Yvectorized(end-IRFw+1:end,:,:) );
        end;
        
        %% preparation of IRF
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
                end
            end
        end
        
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
        sumB      = sum(sum(sum( max(0,ProfileT + ProfileN))));
        
        %% Estimate D
        sumB = sum(sum(sum( max(0,ProfileT + ProfileN))));
        Z_Cleaned  = Z_FextendedReshape;% max(Z_FextendedReshape -  max(0,ProfileT + ProfileN), 0); % N T L R ////////
        clear Z_FextendedReshape
        [vvvvv, Depth] = max(sum(Z_Cleaned,3),[],2); % N 1 1 R
        Depth      = squeeze(Depth); % N R
        
        v =   indGraph(2,:);
        
        %% Estimate Intensities
        if(T<800)
            for r=1:R % multiscale
                %%% Compute intensity
                for n = 1:N
                    %     int = max(1,Depth(n)+minh):min(T,Depth(n)+maxh);
                    %         int        = max(1,Depth(n,r) - 7):min(T,Depth(n,r)+40);
                    for ell=1:L % Wavelength
                        int             = max(1,Depth(n,r) -trailing(ell)):min(T,Depth(n,r)+Attack(ell));
                        weights         = Z_Cleaned(n,int,ell,r);%Z_F(int,n,:,r);%    max(0,Zdump0(int,n) - bb(:) -  sqrt(ddd)*1*sqrt(bb(:)));
                        Intens(n,ell,r) = sum(weights,2);%sum(weights,1);
                    end
                end
            end
            Rscales = reshape(Intens,[row,col,L,R]);
        else
            Rscales  = reshape(sum(Z_Cleaned,2),[row,col,L,R]); % r c  L R
        end
        
        Intens     = reshape(Rscales(:,:,1),row*col,1);
        Depth      = reshape(Depth(:,1),row*col,1);
        
        B(:,:,1)   = reshape(ProfileN,row,col);
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Hist_Back' %% Input shape of Back, use it to estimate D,R
        %%%Var2 = Back
        %%% Estimate D, R
        
        Back  = Var2 ; % T x 1
        [row,col,T,L] = size(Y);
        N = row*col;
        R = 1; % length(I_resol);
        
        %% Expand Y
        Yvectorized    = permute(reshape(Y, [row*col,T,L]), [2,1,3]); % T N L
        IRFw = 0;
        Yvectorized = Yvectorized./(Back+eps); % T N L
        
        %% preparation of IRF
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
        Z_F = zeros(T+2*IRFw,N,L);
        for ell=1:L
            Hf     = fft(h(:,ell)); % T  1
            Z_F(:,:,ell) = ifft( fft(Yvectorized(:,:,ell)) .*Hf  ,'symmetric'); % T x N x L
        end
        clear Yf Yvectorized Hf Y
        % % %         Z_F         = permute(Z_F,[2,1,3]); % N T L
        
        %% Estimate D
        sumB = 0;
        Z_Cleaned  = Z_F;% N T L R ////////
        clear Z_FextendedReshape
        % % %         [vvvvv, Depth] = max(sum(Z_F,3),[],2); % N 1 1 R
        [vvvvv, Depth] = max(sum(Z_F,3)); % N 1 1 R
        Depth    = permute(Depth,[2,3,1]); % N R
        
        B(:,:,1)   = zeros(row,col);
        
        %% Estimate Intensities
        if(T<800)
            for r=1:R % multiscale
                %%% Compute intensity
                for n = 1:N
                    %     int = max(1,Depth(n)+minh):min(T,Depth(n)+maxh);
                    %         int        = max(1,Depth(n,r) - 7):min(T,Depth(n,r)+40);
                    for ell=1:L % Wavelength
                        int             = max(1,Depth(n,r) -trailing(ell)):min(T,Depth(n,r)+Attack(ell));
                        % %                         weights         = Z_F(n,int,ell,r);%Z_F(int,n,:,r);%    max(0,Zdump0(int,n) - bb(:) -  sqrt(ddd)*1*sqrt(bb(:)));
                        weights         = Z_F(int,n,ell,r);%Z_F(int,n,:,r);%    max(0,Zdump0(int,n) - bb(:) -  sqrt(ddd)*1*sqrt(bb(:)));
                        Intens(n,ell,r) = sum(weights,1);%sum(weights,1);
                    end
                end
            end
            Rscales = reshape(Intens,[row,col,L,R]);
        else
            % %             Rscales  = reshape(sum(Z_F,2),[row,col,L,R]); % r c  L R
            Rscales  = reshape(sum(Z_F,1),[row,col,L,R]); % r c  L R
        end
        
        Intens     = reshape(Rscales(:,:,:,1),row*col,1);% N  L 1
        Depth      = reshape(Depth(:,1),row*col,1);
        sumB       = 0;
        
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Estimates'
        [row,col,T,L] = size(Y); T=1;L=1;
        N = row*col;
        R = 1;%
        Depth      = Y(:) ;%reshape(Depth(:,1),row*col,1);
        Intens     = Var2(:);%reshape(Rscales(:,:,1),row*col,1);
        sumB       = 0;
        B(:,:,1)   = zeros(row,col);
        
end
v =   indGraph(2,:);

%%%% output Depth, Intens
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Expand D, R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Int(:,:,:) = reshape(Intens(:,1),[row,col,L])  ;
Dep(:,:,1) =  reshape(Depth(:,1),row,col);
R          =   length(I_resol);

% Intens(Depth(:,1)==0)  =  min(mean(Intens(Depth(:,1)>0) )/2, mean(Intens(Depth(:,1)==0) ))  ; % AH-mod
%Intens(Depth(:,1)==0)  =  min(mean(Intens(Depth(:,1)>0) ), mean(Intens(Depth(:,1)==0) ))  ;  % AH-mod
Intens(Depth(:,1)==0) = 0; 
Int(:,:,:) = reshape(Intens(:,1),[row,col,L])  ; % AH-mod

if(usewmf==1)
    % Weighted Median
    for dim    = 2:R
        weight  = Intens(Neighbours(:,1:v(dim)))./sum(Intens(Neighbours(:,1:v(dim))),2);
        for mm  = [1:-1:1]
            weightD = 1;%abs(Depth(Neighbours(:,1:v(dim)))-ddd(:)) <15*mm;
            weight  = weight.*weightD ./ sum( weight.*weightD,2);
        end
        Dep(:,:,dim)    = WeightedMedian_Parallel(Dep(:,:,1), [], N,weight,Neighbours(:,1:v(dim)), indGraph(:,1:(dim)),1);
        Depth2          = Dep(:,:,dim);
        weightD         = exp(-abs(Depth(Neighbours(:,1:v(dim)))-Depth2(:)))./sum(exp(-abs(Depth(Neighbours(:,1:v(dim)))-Depth2(:)) ) ,2);
        for ell=1:L
            Int(:,:,ell,dim)    = WeightedMedian_Parallel(Int(:,:,ell,1), [], N,weightD,Neighbours(:,1:v(dim)), indGraph(:,1:(dim)),0);
        end
        B(:,:,dim)      = WeightedMedian_Parallel(B(:,:,1), [], N,weightD,Neighbours(:,1:v(dim)), indGraph(:,1:(dim)),1);
    end
    %figure;for i=1:3, subplot(2,3,i);imagesc(reshape(Dep(:,:,i),row,col));subplot(2,3,i+3);imagesc(reshape(Int(:,:,1,i),row,col));end
    %sgtitle('After WMF');
else
    % ToF approximation
    T = max(Depth(:));
    h= normpdf(1:T,T/2+1,15); %%% To tune
    for dim = 2:R
        HistN = zeros(T,N);
        d_neighb = Depth(Neighbours(:,1:v(dim)));
        r_neighb = Intens(Neighbours(:,1:v(dim)));
        r_neighb(d_neighb==0) = 0;
        d_neighb(d_neighb==0) = 1;
        for nnd = 1:N
            HistN(d_neighb(nnd,:),nnd) =  r_neighb(nnd,:);
        end
        HistN2 = conv2(HistN,h(:),'same');
        [vvv, DhrNew]   = max(HistN2);
        Dep(:,:,dim)    = reshape(DhrNew, row,col);
        for nnd = 1:N
            Int2(nnd,dim) = sum(HistN( max(DhrNew(nnd)-10,1):min(DhrNew(nnd)-10,T)  ,nnd));
        end
        Int(:,:,dim) = reshape(Int2(:,dim), row,col);
    end
    %figure;for i=1:3, subplot(2,3,i);imagesc(reshape(Dep(:,:,i),row,col));caxis([700 900]);subplot(2,3,i+3);imagesc(reshape(Int(:,:,i),row,col));end
end

% keyboard 
Dscales   = Dep* Tbin*3*10^8/2;
Rscales    = Int; % r c  L R
D_uncert=zeros(row,col);
R_uncert=zeros(row,col);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Dscales=Dscales(:,:,R);
%Rscales=Rscales(:,:,R);
InMat      = reshape(Rscales,[N  R]); % N x R  
Int      = squeeze(mean(InMat,2)); % N x Lambda

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filter D, R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(useguide==1)
    ParamPC = mean(mean(sum(InMat,2)))>10;
    DnMat      = reshape(Dscales(:,:,1:R),[N R]); % N x R
      
    [ Dguide, Lvect , D_up_bar]   = Build_Dguide_Graph_neighbors_v2(DnMat,row,col,Neighbours,local_shifts,ThreshDep,GuideD,ParamPC,Tbin,CvPC);
    dguide = reshape(Dguide',row,col,[]);
    %figure;for i=1:3, subplot(1,3,i);imagesc(dguide(:,:,i));caxis([700 900]*Tbin*3*10^8/2);end
    % [ I_up_bar, Iguide, Lvect ]   = Build_Iguide_Graph_neighbors(InMat,row,col,Neighbours,GuideI);

    D = zeros(1,N);
    W = zeros(R,N);
    for r=1:R
        ind=find(Dguide(r,:)>0 & D==0);
        D(ind) = Dguide(r,ind);
        W(r,ind) = 1;
    end
    W = (D_up_bar == D);
    sW = sum(W);
    ind = find(sW==0);
    W   = W./sW;
    W(:,ind) = 0;
    Dscales = reshape(D,row,col);
elseif(useguide==2)
    Dscales=median(Dscales,3);
elseif(useguide==3)
    Dscales=mean(Dscales,3);
elseif(useguide==4)
    Dscales=mean(Dscales(:,:,2:end),3);
end

Int     = reshape(Int(:,:,1),row,col);
Rscales = Int.*(Dscales>0);
Dep     = Dep(:,:,1)* Tbin*3*10^8/2;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% uncertainty on D
D_up_bar3 = kron(D_up_bar([1:Lvect:end],:), ones(Lvect,1) ); % nd SCale  x N  noisy D shifted

p = 1;D_up_bar2 = [];Dguide = Dguide';
            for ell=1:R
                for k=1:Lvect
                    D_up_bar2(p,:) =   Dguide(Neighbours(:,k),ell)  ; % nd SCale  X  N    clean D shifted
                    p = p+1;
                end
            end
% Duncert = sum( (D_up_bar2>0)'.*(D_up_bar>0)'.*exp(- abs(D_up_bar - D_up_bar2) /( 2*ThreshDep))'.* abs(D(:) - D_up_bar2'),2) / (1+ size(D_up_bar,1));

Weight    = exp(- abs(D_up_bar3 - D_up_bar2) /( 2*ThreshDep))';
Weight    = Weight./(sum(Weight,2)+eps);
D_uncert  = sum( Weight.* abs(D(:) - D_up_bar'),2) / (1+ size(D_up_bar,1));


%% uncertainty on R 
Int        = squeeze(mean(InMat,3)).*(Dscales(:)>0); % N x L
for ell=1:L
    I   = squeeze(InMat(:,ell,:)).*(Dscales(:)>0); % N   R
    R_uncert(:,ell) = sqrt(mean((Int(:,ell) - I).^2,2)); 
end
%figure;subplot(2,2,2);imagesc(reshape(D_uncert,row,col));colorbar
%subplot(2,2,1);imagesc(reshape(D,row,col));colorbar
%subplot(2,2,3);imagesc(reshape(Rscales,row,col));colorbar
%subplot(2,2,4);imagesc(reshape(R_uncert,row,col));colorbar


