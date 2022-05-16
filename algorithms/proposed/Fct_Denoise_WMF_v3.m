function  [Dscales, Rscales, sumB, D_uncert, R_uncert, W, Dep, Int] = ...
    Fct_Denoise_WMF_v3(Y, Var2, Neighb, pars, F, Dtype,GuideD,GuideI,CvPC, usefilter, useguide)

switch Dtype % check dta format
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Histogram'
        %% input: I_resol, Y, F, IRFw,
        %% output: Dscales, Rscales, Back
        
        [row,col,T,L] = size(Y);
        N = row*col;
        R = length(Neighb.I_resol);
        
        %% Expand Y
        Yvectorized    = permute(reshape(Y, [row*col,T,L]), [2,1,3]); % T N L
        if(T > 800 )
            pars.IRFw = 0;
        else
            Yvectorized   = cat(1,Yvectorized(pars.IRFw:-1:1,:,:),Yvectorized,Yvectorized(end-pars.IRFw+1:end,:,:) );
        end;
        
        %% preparation of IRF
        if(prod(size(F))> (T*L))
            h = squeeze(F(:,floor(T/2),:));
        else
            h = F;
        end
        
        if size(h,1)>(T+2*pars.IRFw)
            h = h(1:sY3,:);
        else
            d = zeros((T+2*pars.IRFw),size(h,2));
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
        Z_F = zeros(T+2*pars.IRFw,N,L);
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
                if(Neighb.I_resol(r)>1)
                    %             Z_Fextended(:,:,:,ell,r)  = imboxfilt3(Z_Fextended(:,:,:,ell,1) , [Neighb.I_resol(r) Neighb.I_resol(r) 1],'Padding' ,'replicate');
                    Z_Fextended(:,:,:,ell,r)   = convn(Z_Fextended(:,:,:,ell,1) , ones(Neighb.I_resol(r),Neighb.I_resol(r))/Neighb.I_resol(r)^2,'same');
                end
            end
        end
        
        if(pars.IRFw>0)
            Z_Fextended       = Z_Fextended(:,:,pars.IRFw+1:T+pars.IRFw,:,:); % r c T  L R
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
        
        v =   Neighb.indGraph(2,:);
        
        %% Estimate Intensities
        if(T<800)
            for r=1:R % multiscale
                %%% Compute intensity
                for n = 1:N
                    %     int = max(1,Depth(n)+minh):min(T,Depth(n)+maxh);
                    %         int        = max(1,Depth(n,r) - 7):min(T,Depth(n,r)+40);
                    for ell=1:L % Wavelength
                        int             = max(1,Depth(n,r) -pars.trailing(ell)):min(T,Depth(n,r)+pars.Attack(ell));
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
        R = 1; % length(Neighb.I_resol);
        
        %% Expand Y
        Yvectorized    = permute(reshape(Y, [row*col,T,L]), [2,1,3]); % T N L
        pars.IRFw = 0;
        Yvectorized = Yvectorized./(Back+eps); % T N L
        
        %% preparation of IRF
        if(prod(size(F))> (T*L))
            h = squeeze(F(:,floor(T/2),:));
        else
            h = F;
        end
        
        if size(h,1)>(T+2*pars.IRFw)
            h = h(1:sY3,:);
        else
            d = zeros((T+2*pars.IRFw),size(h,2));
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
        Z_F = zeros(T+2*pars.IRFw,N,L);
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
                        int             = max(1,Depth(n,r) -pars.trailing(ell)):min(T,Depth(n,r)+pars.Attack(ell));
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
v =   Neighb.indGraph(2,:);

%%%% output Depth, Intens
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Expand D, R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Intens(Depth(:,1)==1)  =  min(mean(Intens(Depth(:,1)>1) ), mean(Intens(Depth(:,1)==1) ))  ; 
Intens(Depth(:,1)<=1) = 0; 
Int(:,:,:) = reshape(Intens(:,1),[row,col,L])  ; % AH-mod
Dep(:,:,1) =  reshape(Depth(:,1),row,col);
R          =   length(Neighb.I_resol);

%run mres_denoiser.m
%sum(sum(Dep(:,:,1)~=dep_preproc))

if(usefilter==1)
    % Weighted Median
    for dim    = 2:R
        weight  = Intens(Neighb.neighb(:,1:v(dim)))./sum(Intens(Neighb.neighb(:,1:v(dim))),2);
        for mm  = [1:-1:1]
            weightD = 1;%abs(Depth(Neighbours(:,1:v(dim)))-ddd(:)) <15*mm;
            weight  = weight.*weightD ./ sum( weight.*weightD,2);
        end
        Dep(:,:,dim)    = WeightedMedian_Parallel(Dep(:,:,1), [], N,weight,Neighb.neighb(:,1:v(dim)), Neighb.indGraph(:,1:(dim)),1);
        Depth2          = Dep(:,:,dim);
        weightD         = exp(-abs(Depth(Neighb.neighb(:,1:v(dim)))-Depth2(:)))./sum(exp(-abs(Depth(Neighb.neighb(:,1:v(dim)))-Depth2(:)) ) ,2);
        for ell=1:L
            Int(:,:,ell,dim)    = WeightedMedian_Parallel(Int(:,:,ell,1), [], N,weightD,Neighb.neighb(:,1:v(dim)), Neighb.indGraph(:,1:(dim)),0);
        end
        B(:,:,dim)      = WeightedMedian_Parallel(B(:,:,1), [], N,weightD,Neighb.neighb(:,1:v(dim)), Neighb.indGraph(:,1:(dim)),1);
    end
    %figure;for i=1:3, subplot(2,3,i);imagesc(reshape(Dep(:,:,i),row,col));subplot(2,3,i+3);imagesc(reshape(Int(:,:,1,i),row,col));end
    %sgtitle('After WMF');
elseif(usefilter==2)
    % Histogram approximation
    T = max(Depth(:));
    h= normpdf(1:T,T/2+1,15); %%% To tune
    for dim = 2:R
        HistN = zeros(T,N);
        d_neighb = Depth(Neighb.neighb(:,1:v(dim)));
        r_neighb = Intens(Neighb.neighb(:,1:v(dim)));
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
elseif(usefilter==3) % ToF approx
    for dim = 2:R    
        %dn=[680, 680, 676,   0, 675, 679, 701,   0,   0, 687, 677,   0,   0,   0, 678,   0, 680, 679,   0,   0, 476, 361,   0, 534, 534];
        %rn=[2,  2,  3,  0,  4,  3,  1,  0,  0,  3,  2,  0,  0,  0,  3,  0,  1,  3,  0,  0,  1,  1,  0, 1, 1];
        %%rn(dn==1)  =  min(mean(rn(dn>1)), mean(rn(dn==1)))  ; 
        d_neighb = Depth(Neighb.neighb(:,1:v(dim)));    
        r_neighb = Intens(Neighb.neighb(:,1:v(dim)));
        %isequal(sort(d_neighb,2), sort(dep_neighb(1:v(dim),:))')
        %r_neighb(d_neighb==1) = 0;    
        %d_neighb(d_neighb==0) = 1;    
        for nnd = 1:N        
            dn = d_neighb(nnd,:);       
            rn = r_neighb(nnd,:)  ;     
            mask = rn==0;
            dn(mask) = [];
            rn(mask) = [];        
            if(isempty(dn))            
                DhrNew3(nnd, dim) = 1;            
                IntTemp(nnd,dim) =0;        
            else           
                DiffV = abs(dn - dn(:));            
                vvvb =  (rn(:)*rn) .* exp( - DiffV.^2/(2*pars.sigIRF^2));            
                [m, indD] = max(sum(vvvb));     %********            
                %DhrNew3(nnd) = dn(indD);     %********            
                IntTemp(nnd,dim) = sum((DiffV(:,indD)<(pars.sigIRF*3)) .*rn(:));       
                DhrNew3(nnd,dim) = sum((DiffV(:,indD)<(pars.sigIRF*3)) .*dn(:).*rn(:))/ (eps+sum((DiffV(:,indD)<(pars.sigIRF*3)) .*rn(:)));
                %IntTemp(nnd,dim)
                %DhrNew3(nnd,dim)
            end    
        end    
        Dep(:,:,dim) = reshape(DhrNew3(:,dim), row,col);    
        Int(:,:,dim) = reshape(IntTemp(:,dim), row,col)/v(dim);    
    end    
end

%a=Dep(:,:,2);
%b=abs(a-dep_filter3);
%sum(abs(a(:)-dep_filter3(:))>0.01)
%[c1, c2]= find(abs(a-dep_filter3)>0.1)
%Dep(:,:,2)=dep_filter3;

%a=Dep(:,:,3);
%b=abs(a-dep_filter5);
%sum(abs(a(:)-dep_filter5(:))>0.1)

% keyboard 
Dscales   = Dep* pars.Tbin*3*10^8/2;
Rscales    = Int; % r c  L R
D_uncert=zeros(row,col);
R_uncert=zeros(row,col);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Dscales=Dscales(:,:,R);
%Rscales=Rscales(:,:,R);
InMat      = reshape(Rscales,[N  R]); % N x R  
Int      = squeeze(mean(InMat,2)); % N x Lambda

%run mres_dguide.m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filter D, R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(useguide==1)
    ParamPC = mean(mean(sum(InMat,2)))>10;
    DnMat      = reshape(Dscales(:,:,1:R),[N R]); % N x R
      
    [ Dguide, Lvect , D_up_bar]   = Build_Dguide_Graph_neighbors_v2(DnMat,row,col,Neighb.neighb,Neighb.local_shifts,pars.ThreshDep,GuideD,ParamPC,pars.Tbin,CvPC);
    dguide = reshape(Dguide',row,col,[]);
    figure;for i=1:3, subplot(1,3,i);imagesc(dguide(:,:,i));caxis([2 15]*pars.Tbin*3*10^8/2);end
    keyboard
    % [ I_up_bar, Iguide, Lvect ]   = Build_Iguide_Graph_neighbors(InMat,row,col,Neighb.neighb,GuideI);
    %sum(sum(abs(dguide(:,:,1)-dguide_res_1)>0.1))
    %sum(sum(abs(dguide(:,:,2)-dguide_res_3)>0.1))
    %sum(sum(abs(dguide(:,:,3)-dguide_res_7)>0.1))
    %a=dguide(:,:,1);
    %b=find(abs(a-dguide_res_1)>0.1);
    
    %[dv,di]=sort(abs(a(:)-dguide_res_1(:)),'descend');

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
Dep     = Dep(:,:,1)* pars.Tbin*3*10^8/2;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% uncertainty on D
D_up_bar3 = kron(D_up_bar([1:Lvect:end],:), ones(Lvect,1) ); % nd SCale  x N  noisy D shifted

p = 1;D_up_bar2 = [];Dguide = Dguide';
            for ell=1:R
                for k=1:Lvect
                    D_up_bar2(p,:) =   Dguide(Neighb.neighb(:,k),ell)  ; % nd SCale  X  N    clean D shifted
                    p = p+1;
                end
            end
% Duncert = sum( (D_up_bar2>0)'.*(D_up_bar>0)'.*exp(- abs(D_up_bar - D_up_bar2) /( 2*pars.ThreshDep))'.* abs(D(:) - D_up_bar2'),2) / (1+ size(D_up_bar,1));

Weight    = exp(- abs(D_up_bar3 - D_up_bar2) /( 2*pars.ThreshDep))';
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


