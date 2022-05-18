function [WR_scale,WRsym_scale] = Compute_Weights_Intens(Iguide,I_up_bar,Lvect,I_resol,r, c,GuideI);
 
[R0,N,L]   = size(Iguide);
R0         = R0/Lvect;
%%%
% I_resol = 1./[1 3 9];%ones(size(I_resol));
% ThreshI = 0.1 + 5* squeeze(Iguide((R0-1)*Lvect,:,:))./max( squeeze(Iguide((R0-1)*Lvect,:,:)));
ThreshI = 0.5*(max(0.1 , squeeze(Iguide((R0-1)*Lvect,:,:)) ));

% ThreshI = 0.1  + 0.5* squeeze(Iguide((R0-1)*Lvect,:,:))./max( squeeze(Iguide((R0-1)*Lvect,:,:)));
  
for Wave = 1:L
    ell=1; 
    if(GuideI==0) % Avoid weights to 1 when not using guide
        I_up_bar(1,:,Wave) = I_up_bar(ell*Lvect+1,:,Wave);
    end
    
    WR     = exp(- ( abs(squeeze(Iguide(ell,:,Wave)    )  ...
        - I_up_bar(Lvect*(ell-1)+1:ell*Lvect,:,Wave) )./(2*ThreshI(:,Wave)'*I_resol(min(R0,ell)))).^2);% ... %;% ...
    WRsym  = exp(- ( abs(squeeze(I_up_bar(ell,:,Wave)    )  ...
        - Iguide(Lvect*(ell-1)+1:ell*Lvect,:,Wave) )./(2*ThreshI(:,Wave)'*I_resol(min(R0,ell)))).^2);% ... %;% ...
    
% %     if(GuideI==0) % Avoid weights to 1 when not using guide 
% %         WR(1,:) =exp(- ( abs(squeeze(Iguide(1,:,Wave) )  ...
% %             - I_up_bar(ell*Lvect+1,:,Wave) )./(2*ThreshI(:,Wave)'*I_resol(min(R0,ell)))).^2);% ... %;% . 
% %     end
    coeffPast  = 1;
    for ell=2:R0-1
        coeffPast = 1;%coeffPast * prod(1-Wscale((ell-2)*Lvect+1,:),1).^(PowerPast);
        %      %%% Compare pixel with neighbours and with median

        WR   = [WR; coeffPast ...
            .*exp(- ( abs(squeeze( Iguide(ell,:,Wave))  - I_up_bar(Lvect*(ell-1)+1:ell*Lvect,:,Wave) )./(2*ThreshI(:,Wave)'*I_resol(min(R0,ell)))).^2 )   ]; %...
        WRsym   = [WRsym; coeffPast ...
            .*exp(- ( abs(squeeze( I_up_bar(ell,:,Wave))  - Iguide(Lvect*(ell-1)+1:ell*Lvect,:,Wave) )./(2*ThreshI(:,Wave)'*I_resol(min(R0,ell)))).^2 )   ]; %...
         
    end
    WR    = [WR; coeffPast ...
        .*exp(-  (abs(squeeze(Iguide(R0,:,Wave))  - I_up_bar(Lvect*(R0-1)+1:R0*Lvect,:,Wave) )./(2*ThreshI(:,Wave)'*I_resol(min(R0,ell)))).^2) ];%L N
    WRsym    = [WRsym; coeffPast ...
        .*exp(-  (abs(squeeze(I_up_bar(R0,:,Wave))  - Iguide(Lvect*(R0-1)+1:R0*Lvect,:,Wave) )./(2*ThreshI(:,Wave)'*I_resol(min(R0,ell)))).^2) ];%L N
    
    WR_scale(:,:,Wave) = WR; % 9 Scales X N X Lambda
    WRsym_scale(:,:,Wave) = WRsym; % 9 Scales X N X Lambda
    clear WR
end
% WRsym_scale=WR_scale;
% WR_scale    = histeq((WR_scale));
% WRsym_scale = histeq((WRsym_scale));

% figure;imagesc(reshape(WR_scale(10,:,1),r,c))
% keyboard