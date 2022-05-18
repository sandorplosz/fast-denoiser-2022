function [Dep, Int3,Int,Int2, Br,Eps_d,Eps_I]= Robust_Median_Bayesian(Dnoisy,Inoisy,I_resol,Neighbours, NeighboursSR...
    ,weightsTemp,local_shifts_Med, local_shifts_WM,Neighbours_WM, indGraph_WM...
    ,ThreshDep,FramHist,Sig2_IRF,NbrePulse,Tbin,SR,GuideD,GuideI)

%% Check the entered parameters; obtain dimensions
%--------------------------------------------------------------
[r c Rpast] = size(Dnoisy);
R0      = length(I_resol);
N       = r*c;
L       = size(Inoisy,3);

DnoisyPC   = reshape(Dnoisy,[N Rpast]);% N R T
DnMat      = reshape(Dnoisy(:,:,1:R0),[N R0]); % N x R
InMat      = reshape(Inoisy(:,:,:,1:R0),[N L R0]); % N x R
Lvect      = size(Neighbours,2);
Neighbours = Neighbours(1:N,:);


%--------------------------------------------------------------
%% Build Spatial Weights: weights for Depth and Reflectivity
%--------------------------------------------------------------
Generate_Multiscale_Weights
%% %% Optional super-resolution, not used in this version since SR=0
rSR = (r*(2^SR )-(2^SR -1));
cSR = (c*(2^SR )-(2^SR -1)) ;
N = rSR * cSR;

if(SR>0)
    clear InMat;keyboard
    for lll=1:L
        for rrr=1:R0
            InMat(:,:,lll,rrr) = interp2(Inoisy(:,:,lll,rrr),SR,'linear');
        end
    end
    InMat = reshape(InMat,[N L R0]);
    
    
    
    Wscale = reshape(Wscale',[r,c,Lvect*R0]);
    Wsym   = reshape(Wsym',[r,c,Lvect*R0]);
    Dguide = reshape(Dguide',[r,c,Lvect*R0]);
    for klk=1:Lvect*R0
        Wscale00(:,:,klk) =   interp2(Wscale(:,:,klk),SR,'linear')/(SR+1)^2  ;   % 9 Scales X N X Lambda
        Wsym00(:,:,klk) =   interp2(Wsym(:,:,klk),SR,'linear')/(SR+1)^2  ;   % 9 Scales X N X Lambda
        Dguide00(:,:,klk) =   interp2(Dguide(:,:,klk),SR,'nearest')  ;   % 9 Scales X N X Lambda
    end
    Wscale = reshape(Wscale00,N,Lvect*R0)' ; weights0=Wscale';
    Wsym   = reshape(Wsym00,N,Lvect*R0)' ; %weights0=Wscale';
    Dguide = reshape(Dguide00,N,Lvect*R0)' ; 
    
    
    % Wscale = reshape( interp2(reshape(Wscale',[r,c,Lvect*R0]),SR,'linear')/SR^2,N,Lvect*R0)'  ;   % 9 Scales X N X Lambda
    I_up_bar0 = reshape(permute(I_up_bar,[2,1,3]),[r,c,Lvect*R0,L]);
    for klk=1:Lvect*R0
        for lll=1:L
            I_up_bar1(:,:,klk,lll) =   interp2(I_up_bar0(:,:,klk,lll),SR,'linear')  ;   % 9 Scales X N X Lambda
        end
    end
    I_up_bar = permute(reshape(I_up_bar1, [N,Lvect*R0,L] ), [2,1,3] ) ;
    % I_up_bar = permute(reshape(interp2(   reshape(permute(I_up_bar,[2,1,3]),[r,c,Lvect*R0,L])    ,SR,'linear'), [N,Lvect*R0,L] ), [2,1,3] );
    
    
    DnoisyPC0 = reshape(DnoisyPC,[r,c,Rpast]);
    for rrr=1:Rpast
        DnoisyPC1(:,:,rrr) =   interp2(DnoisyPC0(:,:,rrr),SR,'linear')  ;   % 9 Scales X N X Lambda
    end
    DnoisyPC = reshape(DnoisyPC1,N,Rpast) ;
    
    % DnoisyPC =  reshape( interp2(reshape(DnoisyPC,[r,c,Rpast]),SR,'nearest'),N,Rpast)  ;
    WR_scale0 = reshape(permute(WR_scale,[2,1,3]),[r,c,Lvect*R0,L]);
    WRsym_scale0 = reshape(permute(WRsym_scale,[2,1,3]),[r,c,Lvect*R0,L]);
    for klk=1:Lvect*R0
        for lll=1:L
            WR_scale1(:,:,klk,lll) =   interp2(WR_scale0(:,:,klk,lll),SR,'linear')/(SR+1)^2  ;   % 9 Scales X N X Lambda
            WRsym_scale1(:,:,klk,lll) =   interp2(WRsym_scale0(:,:,klk,lll),SR,'linear')/(SR+1)^2  ;   % 9 Scales X N X Lambda
        end
    end
    WR_scale = permute(reshape(WR_scale1,[N,Lvect*R0,L]),[2,1,3]) ;
    WRsym_scale = permute(reshape(WRsym_scale1,[N,Lvect*R0,L]),[2,1,3]) ;
    
    rowSR = (r*(2^SR )-(2^SR -1));
    colSR = (c*(2^SR )-(2^SR -1));
%     [Neighbours_MedSR indGraph_MedSR,local_shifts_MedSR] =  Build_Graph_Neighbours_Array_v2(rowSR,colSR,[1 floor(WindSR/2)*((2^SR -1)+1)*2+1]);% Define graph of correlations between pixels
  
% % %     [NeighboursSR indGraph_SR,local_shifts_SR] =  Build_Graph_Neighbours_Array_v2(rowSR,colSR,[1 3]);% Define graph of correlations between pixels
    [NeighboursSR indGraph_SR,local_shifts_SR] =  Build_Graph_Neighbours_Array_v2(rowSR,colSR,[1 5]);% Define graph of correlations between pixels
    NeighboursSR  = NeighboursSR(1:end-1,:);
    N = size(NeighboursSR,1);
    Wscale2      = permute(reshape(Wsym,[Lvect,R0,N]),[1,3,2]); % nd N R  //////
    
end




%--------------------------------------------------------------


%% Reflectivity estimate
%--------------------------------------------------------------
% keyboard
PriorI     = squeeze(mean(mean(var(Inoisy(:,:,:,:),[],4))));
aPrior     = 2;
bPrior     = (aPrior+1)*PriorI/10;
bPrior     = permute(repmat(bPrior,[1, Lvect*R0,N]),[2,3,1]);
IterMax    = 4; % maximum number of iterations
Int2       = squeeze(mean(InMat,3)); % N x Lambda
Int        = squeeze(sum(WR_scale.*I_up_bar,1)); %3);
I_up_bar2  = I_up_bar;
% sumWR      = sum(WR_scale,1);
% WR_scale   = WR_scale./sum(WR_scale,1);
DenomEps   = (1+size(WRsym_scale,1)/2 + aPrior);
   
ThconvI    = 0.001;

i = 1;
while(i <= IterMax)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%      Step 1: Update M     %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     M        = (sum(WR_scale.*I_up_bar2,1)./(sumWR));  % 1 x N x Wave GaussMRF
    M        = (sum(WR_scale.*I_up_bar2,1));  % 1 x N x Wave GaussMRF
%     figure(10);subplot(2,2,min(i,4));imagesc(reshape(M,[r,c,3])/max(M(:)));title('Proposed Intensity')
    
    Mtot(:,:,i) = M;
    for  Wave = 1:L
        vtempp            = squeeze(M(1,:,Wave));
        Mextend(:,:,Wave) = vtempp(NeighboursSR)';% nd N Wave
    end
    Mextend=max(eps,Mextend); % nd N Wave = 9 N 3
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%  Step 2: Update  Eps_I    %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    diffI  =  WRsym_scale.*(M - I_up_bar2).^2/2  + bPrior;   % scale nd x N x Wave
%     Eps_I  = squeeze(max(eps,sum(diffI  ,1,'omitnan' ) ./(1+size(WRsym_scale,1)/2 + aPrior)));  % N x Wave
    Eps_I  =  (max(eps,sum(diffI  ,1,'omitnan' ) ./DenomEps));  % 1 x N x Wave
    %     figure(11);for jj=1:3 subplot(1,3,jj);imagesc(reshape(Eps_I(:,jj),[r,c]));title('STD I');end
%         keyboard
    for  Wave = 1:L
        vtempp = squeeze(Eps_I(1,:,Wave));vtempp=vtempp(:);
        Eps_Iextend(:,:,Wave) = vtempp(NeighboursSR)';% nd N Wave
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Step 3: Update I_up_bar2  %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for rr = 1:R0
        Sig_r =  1./  max(eps,sum(WRsym_scale((rr-1)*Lvect+1:rr*Lvect,:,:) ./Eps_Iextend,1));    % 1 N  Wave
        Mu_r  =  sum(Mextend.* WRsym_scale((rr-1)*Lvect+1:rr*Lvect,:,:) ./Eps_Iextend, 1) .*Sig_r ; % 1 N  Wave
        Mu_Sig_r     = Mu_r - Sig_r;
%         Int3(rr,:,:) = (Mu_Sig_r + sqrt( Mu_Sig_r.^2 + 4*Sig_r.* I_up_bar((rr-1)*Lvect+1,:,:) ) )/2; % 1 N  Wave
        Int3  = (Mu_Sig_r + sqrt( Mu_Sig_r.^2 + 4*Sig_r.* I_up_bar((rr-1)*Lvect+1,:,:) ) )/2; % 1 N  Wave
         
        for Wave=1:L
%             vtempp =  Int3(rr,:,Wave);
            vtempp =  squeeze(Int3(1,:,Wave));vtempp=vtempp(:);
            I_up_bar2((rr-1)*Lvect+1:rr*Lvect,:,Wave)  = vtempp(NeighboursSR)'   ;% nd N Wave
        end
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%     Check convergence     %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(i>1)
        Mcost(i-1)   = sum(sum(abs(Mtot(:,:,i)-Mtot(:,:,i-1))))/sum(sum(abs(Mtot(:,:,i-1))));
        if(Mcost(i-1) < ThconvI)   i = IterMax;end
    end
    i=i+1;
end
Int3 = squeeze(M);
 

%--------------------------------------------------------------
%% Depth estimate using coordinate descent
%--------------------------------------------------------------
SigLike  = Sig2_IRF ./ (InMat/NbrePulse)*( Tbin*3*10^8/2)^2; % % N xL  x R
SigLike  = squeeze(1./(sum(1./SigLike ,2))); % N x R
Eps_d    = 10*mean(SigLike,2) ; %10*mean(SigLike,2) ; % N x 1 0.001
D_MAP    = Dguide(1:Lvect:end,:)'; % DnoisyPC   ;% N R T
DnoisyPC = Dguide(1:Lvect:end,:)'; % change of rough estimates
Thconv   = 0.01; % Threshold to check convergence
% keyboard
i=1;
while(i <= IterMax) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%  Weighted median   %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [x_MAP,  vectMat, weights0] = WeightedMedian_SpaceTime_Bayesian(D_MAP(:,1:R0),D_MAP(:,R0+1:end),weights0,Eps_d,NeighboursSR,NeighboursSR);
    x_MAPTot(:,i) = x_MAP(:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%      Update variances     %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    diffD  = abs(x_MAP(:) - vectMat) ; % N x L
    %     Eps_d  = max(eps,sum(diffD .* weights0  ,2,'omitnan' ) /(1+size(weights0,2)));  % N x 1
    Eps_d  = max(eps,sum(diffD .* Wsym'  ,2,'omitnan' ) ./(1+size(Wsym,1)));  % N x 1
    %     figure(105);imagesc(reshape(Eps_d,r,c));keyboard
    epsTot(:,i) = Eps_d(:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% soft Thresholding  %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %     Run every wavelength independently or
    %         for ell=1:R0
    %           [ D_MAP(:,:,ell)]   = Generalized_Soft_Thresholding(Dnoisy(:,:,ell), SigLike(:,ell), x_MAP, Wscale((ell-1)*Lvect+1:ell*Lvect,:)', Eps_d(Neighbours), 1:N, Neighbours);
    %         end
    % %     Run all wavelengths in parallel
    [D_MAP ]= Generalized_Soft_Thresholding_ParaScales(DnoisyPC(:,1:R0), SigLike , x_MAP(:), Wscale2, Eps_d(NeighboursSR), NeighboursSR);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%     Check convergence     %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(i>1)
        Xcost(i-1)   = sum(abs(x_MAPTot(:,i)-x_MAPTot(:,i-1)))/sum(abs(x_MAPTot(:,i-1)));
        Epscost(i-1) = sum(abs(epsTot(:,i)-epsTot(:,i-1)))/sum(abs(epsTot(:,i-1)));
        if(Xcost < Thconv)   i = IterMax;end
    end
    i = i+1;
end

Dep  = x_MAP(:);
Br   = (reshape(interp2(Dnoisy(:,:,1),SR,'nearest'),[N 1]) - Dep(:))' ; % R x N


return

