function [NormW, W_D_to_X, Wsym, weights0] = Compute_Weights_Depth(D_med2,D_up_bar, ThreshDep, I_resol, R0, Lvect,NeighboursSR,weightsTemp,Neighbours_WM)
%% Input: D_med2 D_up_bar  ThreshDep I_resol R0 Lvect
%% Output:  NormW     W_D_to_X   Wsym     

% compare each pixel of clean data with a window of observed data
PowerPast = mean(mean((abs(D_up_bar(1,:)  - D_med2(1,:)  ) < 2*ThreshDep)  & (D_up_bar(1,:)>0) & (D_med2(1,:) >0) ))>0.5;
ell=1;PowerPast=1;
ValidDmed = D_med2>0;ValidDbar  = D_up_bar>0;
W_D_to_X =ValidDmed(Lvect*(ell-1)+1,:).*ValidDbar(Lvect*(ell-1)+1:ell*Lvect,:).*...
    exp(- ( abs(squeeze(D_med2(Lvect*(ell-1)+1,:))  - D_up_bar(Lvect*(ell-1)+1:ell*Lvect,:) )/(2*ThreshDep*I_resol(min(R0,ell)))).^2 ... %;% ...
    -  (1-PowerPast)*abs(squeeze(D_med2(Lvect*(ell-1)+1,:)    - D_up_bar(Lvect*(ell)+1:(ell+1)*Lvect,:) )/(2*ThreshDep*I_resol(min(R0,ell)))).^2);
Wsym = ValidDbar(Lvect*(ell-1)+1,:).*ValidDmed(Lvect*(ell-1)+1:ell*Lvect,:).* ...
    exp(- ( abs(squeeze(D_up_bar(Lvect*(ell-1)+1,:))  - D_med2(Lvect*(ell-1)+1:ell*Lvect,:) )/(2*ThreshDep*I_resol(min(R0,ell)))).^2);

coeffPast  = 1;
for ell=2:R0-1
    vectTemp  = sort(W_D_to_X((ell-2)*Lvect+1:(ell-1)*Lvect,:),'descend');
    coeffPast = coeffPast .* (1-mean(vectTemp(1:3,:),1)).^(PowerPast);
    
    PowerPast = 1;
    %      %%% Compare pixel with neighbours and with median
    W_D_to_X   = [W_D_to_X; ValidDmed(Lvect*(ell-1)+1,:).*ValidDbar(Lvect*(ell-1)+1:ell*Lvect,:).* coeffPast ...
        .*exp(- ( abs(squeeze(D_med2(Lvect*(ell-1)+1,:))  - D_up_bar(Lvect*(ell-1)+1:ell*Lvect,:) )/(2*ThreshDep*I_resol(min(R0,ell)))).^2  ...
        - (1-PowerPast)* abs(squeeze(D_med2(Lvect*(ell-1)+1,:)    - D_up_bar(Lvect*(ell)+1:(ell+1)*Lvect,:) )/(2*ThreshDep*I_resol(min(R0,ell)))).^2)];
    Wsym   = [Wsym; ValidDbar(Lvect*(ell-1)+1,:).*ValidDmed(Lvect*(ell-1)+1:ell*Lvect,:).* coeffPast(NeighboursSR(:,:))' ...
        .*exp(- ( abs(squeeze(D_up_bar(Lvect*(ell-1)+1,:))  - D_med2(Lvect*(ell-1)+1:ell*Lvect,:) )/(2*ThreshDep*I_resol(min(R0,ell)))).^2 )];
end

vectTemp  = sort(W_D_to_X((R0-2)*Lvect+1:(R0-1)*Lvect,:),'descend');
coeffPast = coeffPast .* (1- mean(vectTemp(1:3,:),1)).^(PowerPast);

%         coeffPast = coeffPast .* prod(1-max(Wscale((R0-2)*Lvect+1:(R0-1)*Lvect,:),[],1),1).^(PowerPast);
W_D_to_X    = [W_D_to_X; ValidDmed(Lvect*(R0-1)+1,:).*ValidDbar(Lvect*(R0-1)+1:R0*Lvect,:).* coeffPast ...
    .*exp(-  (abs(squeeze(D_med2(Lvect*(R0-1)+1,:))  - D_up_bar(Lvect*(R0-1)+1:R0*Lvect,:) )*2/(2*ThreshDep*I_resol(min(R0,ell)))).^2) ];%L N
Wsym    = [Wsym; ValidDbar(Lvect*(R0-1)+1,:).*ValidDmed(Lvect*(R0-1)+1:R0*Lvect,:).* coeffPast(NeighboursSR(:,:))' ...
    .*exp(-  (abs(squeeze(D_up_bar(Lvect*(R0-1)+1,:))  - D_med2(Lvect*(R0-1)+1:R0*Lvect,:) )*2/(2*ThreshDep*I_resol(min(R0,ell)))).^2) ];%L N

%         keyboard
indd = sum(W_D_to_X,1)==0;
W_D_to_X(1:Lvect:end,indd) = 1;
Wsym(1:Lvect:end,indd) = 1; % Sometimes neighbor of n is n but we do not consider this


NormW      = sum(W_D_to_X,1); % N x 1
W_D_to_X   = W_D_to_X./ NormW ; %L N
Wsym       = Wsym./ repmat(NormW(NeighboursSR)',R0,1) ; %L N


nd          = size(Neighbours_WM,2);
weightsTemp = 0.1*weightsTemp(1,[nd+1:end])/sum(weightsTemp(1,[nd+1:end]));

fact          = sum(weightsTemp)+sum(W_D_to_X,1); 
W_D_to_X        = sqrt(2*W_D_to_X./fact) ; %L N
weightsTemp   = sqrt(2*weightsTemp./fact') ; %N x L
weights0      =  [W_D_to_X' , weightsTemp];
Wsym          = sqrt(2*Wsym./sum(Wsym,1)) ; %L N ////////



