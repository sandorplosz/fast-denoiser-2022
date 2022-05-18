function [D_MAP] = Generalized_Soft_Thresholding_ParaScales(D_ML, Sig_IRF_R, x_MAP, Weights, Eps_d, Neighbours_Sup);
%%% D_ML Depth from histograms  N x R
%%% Sig_IRF_R pixelwise variance N x R . Related to IRF width and Counts
%%% x_MAP High-Resolution Depth  N x 1
%%% Weights of size Nx nd R
%%% Eps_d  of size Nhr x Lw
%%% Neighbours_Sup of size Nhr x Lw

[N, R0]         = size(D_ML);
D_MAP           = zeros(N,1);
Lw              = size(Eps_d,2); % Number of neighbours
Sig_IRF_R       = permute(Sig_IRF_R,[3 1 2]) ;% 1 x N x  R0

%%% Assign variables
x0    = permute(D_ML,[1 3 2]); % N x 1 x R0
Mu    = x_MAP(Neighbours_Sup); % N x Lw
alpha = Weights.* Sig_IRF_R ./Eps_d'; %Lw N R
clear D_ML Neighbours_Sup x_MAP Weights Sig_IRF_R Eps_d
%%% Apply weighted median for missing D (i.e., D=0) 
pixInd = (permute(repmat(sum(Mu,2)==0,[1,1,R0]) & x0==0, [2, 1, 3] ));%Lw N R
alpha  = alpha.* (10^8*pixInd);

%%% Compute thresholds
[vectOrderMat IndMat] = sort(Mu,2); % N x Lw   ascend
IndMat2               = IndMat' + repmat((0:Lw:(N-1)*Lw),Lw,1); % Lw  N
WM2                   = alpha(repmat(IndMat2(:),1,R0)+ (0:R0-1)*(N*Lw)   ) ;
alphaSort             = permute(reshape(WM2,[Lw,N,R0]),[2,1,3]) ;% N  Lw R
Mu             = vectOrderMat; % N x Lw
clear IndMat IndMat2 WM2 alpha vectOrderMat
CumWeightsRev  = cumsum(alphaSort,2,'reverse'); % N  Lw  R
CumWeights         = cumsum(alphaSort,2); % N  Lw R
clear alphaSort 
Thresh(:,1,:)      = -CumWeightsRev(:,1,:);
Thresh(:,2*Lw,:)   = CumWeights(:,end,:);
% keyboard
Thresh(:,2:(2*Lw-1),:)   =  repelem(CumWeights(:,1:end-1,:)-CumWeightsRev(:,2:end,:),1,2,1) ; %N x 2 Lw x R
 
Region             = squeeze(sum( repelem(x0-Mu,1,2,1) >Thresh, 2)+1); % N x  R value from 1 to 2 Lw+1
 
Coeff                 = zeros(N,2*Lw+1,R0);
Coeff(:,1,:)          = x0+  CumWeightsRev(:,1,:); %N  Lw R   L+1
Coeff(:,3:2:2*Lw-1,:) = x0+  CumWeightsRev(:,2:Lw,:) - CumWeights(:,1:Lw-1,:); %N  Lw R   L+1
Coeff(:,2*Lw+1,:)     = x0 - CumWeights(:,end,:); %N  Lw R   L+1
Coeff(:,2:2:2*Lw,:)   = repmat(Mu,[1,1,R0]); %N  Lw R     L

%%% Select the values for all pixels
for rr=1:R0
    for i=1:N
        D_MAP(i,rr) = Coeff(i,Region(i,rr),rr); %N  R
    end
end


