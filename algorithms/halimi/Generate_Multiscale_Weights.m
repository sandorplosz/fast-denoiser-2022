%% input: DnMat ,InMat,Neighbours,Lvect,R0,ThreshDep, Neighbours_WM, weightsTemp

%% %%% Build shifted versions of Depths and Intensities
%%% Inputs: DnMat NeighboursSR  Neighbours Lvect R0 N ThreshDep InMat
%%% Outputs: D_up_bar, I_up_bar:shifted versions of depth and Intensity
%%%          Dguide: cleaned depths
% % % % switch GuideD
% % % %     case 0 % Multiscale median
% % % tic
% % %         [ D_up_bar, Dguide, I_up_bar, Lvect ]   = Build_Dref_Iref_Graph_neighbors(DnMat,InMat,NeighboursSR,Neighbours,ThreshDep);
% % % toc
% % % keyboard
% % % tic
% % % %     case 1 % PC denoise
% % %         [ D_up_bar, Dguide, I_up_bar, Lvect ]   = Build_DrefPCdenoise_Iref_Graph_neighbors(DnMat,InMat,r,c,NeighboursSR,Neighbours,ThreshDep);
% % % % end
% % % toc
%  keyboard
ParamPC = mean(mean(sum(InMat,2)))>10;
[ D_up_bar, Dguide, Lvect ]   = Build_Dguide_Graph_neighbors(DnMat,r,c,NeighboursSR,Neighbours,ThreshDep,GuideD,ParamPC,Tbin);
 
[ I_up_bar, Iguide, Lvect ]   = Build_Iguide_Graph_neighbors(InMat,r,c,NeighboursSR,Neighbours,GuideI);
% keyboard

%% Optional SR: 

  
%% Generation of multiscale depth weights
[NormW, Wscale, Wsym, weights0] = Compute_Weights_Depth(Dguide,D_up_bar, ThreshDep, I_resol, R0, Lvect,NeighboursSR,weightsTemp,Neighbours_WM);
 
%% Generation of multiscale intensity weights
[WR_scale,WRsym_scale] = Compute_Weights_Intens(Iguide,I_up_bar,Lvect,I_resol,r, c,GuideI);
 



%% Update R weights using depth weights 
WR_scale0    = WR_scale;
WR_scale     = sqrt(WR_scale)    .* repmat(Wscale,[1,1,L]) ;
WRsym_scale  = sqrt(WRsym_scale) .* repmat(Wsym,[1,1,L]) ;

NormW        = max(eps,sum(WR_scale,1)); % N x 1
WR_scale     = WR_scale./ NormW ; %L N
WRsym_scale  = min(1,WRsym_scale./ repmat(NormW(NeighboursSR)',R0,1) ); %L N

Wscale2      = permute(reshape(Wsym,[Lvect,R0,N]),[1,3,2]); % nd N R  //////



