clear all;  clc;
%%  Demo MATLAB code
%%  Demo MATLAB code
% Code by A. Halimi and S. Plosz
% published in IEEE 2022.
% Contact: A. Halimi,  a.halimi@hw.ac.uk
%   
%%%%% Set constants
BackTitle  = {['Expo. Back'],['Gamma. Back']};
BackTitle2 = {['Expo_Back'],['Gamma_Back']}; 
T       = 2500;% 1000;%
SBR     = logspace(-2,1,10);%logspace(-2,2,30);
Photons = unique(round(logspace(0,3,10)));% unique(round(logspace(0,3,30)));
MC      = 50; 
sigma2 =   40^2;  % 36

Alpha  = 1;
Beta   = 1;
Kf     = 0.001/T; % mean 0.1/T  var=10/T
Theta  = 0.01;

NbrePhoton  = 10; %11;
B  = beta(Alpha,Beta);
B2 = betaln(Alpha,Beta);
limitC   = (nchoosek(NbrePhoton,round(NbrePhoton/2)));


%%%% Coeff online processing + difference of photons
Xph    = 0:NbrePhoton-1;
Coeff = betaln(Alpha+NbrePhoton-Xph,Beta+Xph) - Xph*log(T) ...
    -0.5*log((NbrePhoton-Xph)) - ((NbrePhoton-Xph-1)/2)*log((2*pi*sigma2))-B2-log(T);
Coeff(NbrePhoton+1) = betaln(Alpha,Beta+NbrePhoton)-log(T)*(NbrePhoton-1)-B2-log(T);
for i= 0:NbrePhoton-1
    Possib{i+1} = nchoosek(1:NbrePhoton,NbrePhoton-i);
end
%%%% Coeff Approx Likelihood
Xph    = 0:NbrePhoton-1;
CoeffAppL = betaln(Alpha+NbrePhoton-Xph,Beta+Xph) - Xph*log(T) ...
    -B2-log(T);
CoeffAppL(NbrePhoton+1) = betaln(Alpha,Beta+NbrePhoton)-log(T)*(NbrePhoton-1)-B2-log(T);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%     Generate data  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Create an IRF
h= 1/sqrt(2*pi*sigma2)* exp(- ((1:T)-round(T/2)).^2/2/(sigma2));
h = h(:)/sum(h);
[~,attack] = max(h);
h = circshift(h,-attack);

for Choice_Bshape =  1:2; % 1 Exponential,  2 gamma
    Z_GT = zeros(length(SBR),length(Photons),MC);
    Y    = zeros(length(SBR),length(Photons),T,MC);
    for Choice_Proba  =   1:2; % 1 Proba of true detection,  2 false alarm
        
        for m=1:(MC)
            for i=1:length(SBR)
                for j=1:length(Photons)
                    %                     Ns      = (SBR(i)/(SBR(i)+1)*Photons(j));
                    %                     Prob    = Ns-floor(Ns); u = (rand<Prob);
                    %                     Ns      = floor(Ns)*(u==1) + ceil(Ns)*(u==0);
                    Ns      = ceil((SBR(i)/(SBR(i)+1)*Photons(j)));
                    Nb      = Photons(j) - Ns;
                    dep(i,j,m) =  round(T/4+120  + (T/2-240)*rand(1));
                    Sphoton = max(1,round(dep(i,j,m)+ sqrt(sigma2)*randn(Ns,1)));
                    %                     Bphoton = max(1,round(T*rand(Nb,1)));    Back = 1/T* ones(1,T);% exp(-(1:T)/(T/4)); Back = Back/sum(Back);
                    
                    switch Choice_Bshape
                        case 1
                            v=[];mmm=0;
                            while (mmm<Nb)  vProp=exprnd(T/4 ,[1,1]);
                                if(vProp < T) v=[v; vProp];  end;mmm=length(v);
                            end
                            Back = exp(-(1:T)/(T/4)); Back = Back/sum(Back);
                        case 2
                            MM = 0.4*T;   SS = (0.3*T);  param1 = MM^2/SS^2; param2=SS^2/MM;
                            %                             param1 = 9/4; param2=2*T/9;
                            v=[];mmm=0;
                            while (mmm<Nb)  vProp=gamrnd(param1,param2);
                                if(vProp < T) v=[v; vProp];
                                end;mmm=length(v);
                            end
                            Back = gampdf(1:T,param1,param2);%(1:T).*(param1-1).*exp(-(1:T)/param2) ; Back = Back/sum(Back);
                    end
                    
                    Bphoton = max(1,round(v));
                    tof     = [Sphoton;Bphoton];
                    Im{i,j,m} = tof(randperm(length(tof)));
                    if(Ns>0)
                        Z_GT(i,j,m) = 1;
                    end
                    
                    for k=1:length(tof)
                        Y(i,j,tof(k),m) = Y(i,j,tof(k),m)+1;
                    end
                    
                end
            end
        end
        
        switch Choice_Proba
            case 1
                %%%%%%%%%%%%%%%%%%  PD
                clear Z_Gauss Z_Histo Z_Gauss8 Z_Gauss5
                for m = 1:MC
                    m
                    for i =     1:length(SBR)
                        for j=   1:length(Photons)
                            %                     figure;plot(squeeze(Y(i,j,:,m)))
                            %%%% Unif + Gauss
                            tic
                            %                             [Z_Gauss(i,j,m), p11, p00,p1,p1Spat,logp1Fct, logp0Fct,aprev,C] = detect_AH_v5(Im{i,j,m} ,sigma2,Alpha,Beta,T,NbrePhoton,0.5,limitC,2);
                            [Z_Gauss(i,j,m), DiffLogP(i,j,m), Depth(i,j,m),Coeff22,aprev,p1,p0] = detect_AH_NonUnifBack_v7(Im{i,j,m}   ,sigma2,Alpha,Beta,T,NbrePhoton,0.5,limitC, 3, 1/T*ones(T,1));%5)
                            t_Gauss(i,j,m) = toc;
                            
                            %%%% From C:\Users\ah64\OneDrive - Heriot-Watt University\Teaching\PhD_Amir\Codes\Detection_AH
                            tic
                            [Z_Gauss8(i,j,m), DiffLogP2(i,j,m), Depth2(i,j,m),Coeff2,aprev2,p1,p0] = detect_AH_NonUnifBack_v7(Im{i,j,m}   ,sigma2,Alpha,Beta,T,NbrePhoton,0.5,limitC, 3, Back(:));%5)
                            t_Gauss8(i,j,m) = toc;
                            
                            [Z_Gauss5(i,j,m), p11, p00,p1,p1Spat,logp1Fct, logp0Fct,aprev,C] = detect_AH_v5(Im{i,j,m} ,sigma2,Alpha,Beta,T,5,0.5,(nchoosek(5,round(5/2))),2);
                            t_Gauss5(i,j,m) = toc;
                            
                            
                            %%%% Detect Julian
                            y =  squeeze(Y(i,j,:,m));max_r =  Photons(j);
                            tic
                            %                    [Z_Histo(i,j,m), p11, p00,p1,p1Spat,logp1Fct, logp0Fct,aprev,C] = detect_AH_v3(Im{i,j,m} ,sigma2,Alpha,Beta,T,NbrePhoton,0.5,limitC,0);
                            %                     [Z_Histo(i,j,m), p11, p00,p1,p1Spat,logp1Fct, logp0Fct,aprev,C] = detect_AH_v4(Im{i,j,m} ,sigma2,Alpha,Beta,T,NbrePhoton,0.5,limitC);
                            [Z_Histo(i,j,m),W(i,j)] = detect(y,h,max_r);
                            
                            %                     [Z_Histo( i,j,m),W(i,j)] = detect_JulianTOF(y,h,max_r,tof);
                            %                                         Z_Histo(i,j,m)=1;
                            t_Histo(i,j,m)=toc;
                        end
                        %                 pause
                    end
                end
                Z_Histo  = (Z_Histo>0);
                PD_Gauss = mean(Z_Gauss.*Z_GT,3);
                PD_Gauss8 = mean(Z_Gauss8.*Z_GT,3);
                PD_Gauss5 = mean(Z_Gauss5.*Z_GT,3);
                PD_Histo = mean(Z_Histo.*Z_GT,3);
                
                %%%% PD maps
                figure(1);
                subplot(2,3,1+(Choice_Bshape-1)*3);contourf(log10(Photons),log10(SBR),PD_Histo);
                ylabel({['\fontsize{18}' BackTitle{Choice_Bshape}],'\fontsize{14}SBR'})
                xlabel('Photons','fontsize',14);title('HTD','FontSize',18);caxis([0 1])
                xticks([1 2 3]);xticklabels({'10^1','10^2','10^3'})
                yticks([-2:2]); yticklabels({'10^{-2}','10^{-1}','10^{0}','10^1','10^2'});colorbar; %caxis([-0.2, 0.2])
                subplot(2,3,2+(Choice_Bshape-1)*3);contourf(log10(Photons),log10(SBR),PD_Gauss);  ylabel('SBR','FontSize',14); xlabel('Photons','FontSize',14);title(' ETD','FontSize',18);%caxis([0 1])
                xticks([1 2 3]);xticklabels({'10^1','10^2','10^3'})
                yticks([-2:2]); yticklabels({'10^{-2}','10^{-1}','10^{0}','10^1','10^2'});colorbar; %caxis([-0.2, 0.2])
                subplot(2,3,3+(Choice_Bshape-1)*3);contourf(log10(Photons),log10(SBR),PD_Gauss8);  ylabel('SBR','FontSize',14); xlabel('Photons','FontSize',14);title('GETD','FontSize',18);caxis([0 1])
                xticks([1 2 3]);xticklabels({'10^1','10^2','10^3'})
                yticks([-2:2]); yticklabels({'10^{-2}','10^{-1}','10^{0}','10^1','10^2'});colorbar
                
%                 cd('./results')
%                 save(['Results_Synth_PD_' BackTitle2{Choice_Bshape} '.mat'],'Im','Y','Z_GT','SBR','Photons','Z_Gauss','Z_Histo','t_Gauss','t_Histo','PD_Gauss','PD_Histo','t_Gauss5','PD_Gauss5','t_Gauss8','PD_Gauss8')
%                 saveName = ['Detection_TruePositive_HTD_ETD_GETD_SBR_Photons'];
%                 %                                 saveName = ['Detection_TruePositive_HTD_ETD_GETD_SBR_Photons_Sig6_T1000'];
%                 saveas(gcf, saveName, 'fig');
%                 print(saveName,'-dpng','-r512')
%                 cd ..
                
                % %         %%%% Computational times
                % %         figure;
                % %         subplot(2,1,1)
                % %         loglog(Photons,mean(mean(t_Gauss,1),3)*1000,'b','LineWidth',1.5);hold on
                % %         loglog(Photons,mean(mean(t_Gauss8,1),3)*1000,'g+-','LineWidth',1.5);
                % %         %         loglog(Photons,mean(mean(t_Gauss5,1),3)*1000,'mo-','LineWidth',1.5);
                % %         hold on;loglog(Photons,mean(mean(t_Histo,1),3)*1000,'r','LineWidth',1.5);
                % %         %         legend('ETD (M=10)','ETD (M=8)','ETD (M=5)','HTD')
                % %         legend('ETD (M=10)','ETD (M=8)','HTD')
                % %         ylabel('Time (ms)'); xlabel('Photons')
                % %         set(gca,'FontSize',14);grid
                % % %         cd ./images
                % % %         saveName = ['ComputTime_Prop_Histo_Photons_v2'];
                % % %         saveas(gcf, saveName, 'fig');
                % % %         print(saveName,'-dpng','-r512')
                % % %         cd ..
                
                % %          pause
                % %          pause
            case 2
                %%%%%%%%%%%%%%%%%%  PFA
                clear Z_Gauss Z_Histo Z_Gauss8 Z_Gauss5
                %                 MC    = 50;%
                Z_GT  = zeros(length(Photons),MC);
                Yb    = zeros(length(Photons),T,MC);
                %         Back = exp(-(1:T)/(T/4)); Back = Back/sum(Back);%Back = 1/T*ones(1,T);
                for m=1:(MC)
                    for j=1:length(Photons)
                        % %                         Bphoton = max(1,round(T*rand(Photons(j),1)));
                        switch Choice_Bshape
                            case 1
                                v=[];mmm=0;
                                while (mmm<Nb)  vProp=exprnd(T/4,[1,1]);
                                    if(vProp < T) v=[v; vProp];  end;mmm=length(v);
                                end
                                Back = exp(-(1:T)/(T/4)); Back = Back/sum(Back);
                            case 2
                                MM = 0.4*T;   SS = (0.3*T);  param1 = MM^2/SS^2; param2=SS^2/MM;
                                %                             param1 = 9/4; param2=2*T/9;
                                v=[];mmm=0;
                                while (mmm<Nb)  vProp=gamrnd(param1,param2);
                                    if(vProp < T) v=[v; vProp];
                                    end;mmm=length(v);
                                end
                                Back = gampdf(1:T,param1,param2);%(1:T).*(param1-1).*exp(-(1:T)/param2) ; Back = Back/sum(Back);
                        end
                        Bphoton = max(1,round(v));
                        
                        tof     = [Bphoton];
                        ImB{j,m} = tof(randperm(length(tof)));
                        Z_GT(j,m) = 0;
                        
                        for k=1:length(tof)
                            Yb(j,tof(k),m) = Yb(j,tof(k),m)+1;
                        end
                        
                    end
                end
                
                % %% Processing
                for m = 1:MC
                    m
                    for j=1:length(Photons)
                        %%%% Unif + Gauss
                        tic
                        %                         [Z_Gauss(j,m), p11, p00,p1,p1Spat,logp1Fct, logp0Fct,aprev,C] = detect_AH_v5(ImB{j,m} ,sigma2,Alpha,Beta,T,NbrePhoton,0.5,limitC,3);
                        [Z_Gauss(j,m), DiffLogP(j,m), Depth(j,m),Coeff22,aprev] = detect_AH_NonUnifBack_v7(ImB{j,m}   ,sigma2,Alpha,Beta,T,NbrePhoton,0.5,limitC, 3, 1/T*ones(T,1));%5)
                        t_Gauss(j,m) = toc;
                        
                        tic
                        [Z_Gauss8(j,m), DiffLogP2(j,m), Depth2(j,m),Coeff2,aprev2] = detect_AH_NonUnifBack_v7(ImB{j,m}   ,sigma2,Alpha,Beta,T,NbrePhoton,0.5,limitC, 3, Back(:));%5)
                        t_Gauss8(j,m) = toc;
                        
                        tic
                        [Z_Gauss5(j,m), p11, p00,p1,p1Spat,logp1Fct, logp0Fct,aprev,C] = detect_AH_v5(ImB{j,m} ,sigma2,Alpha,Beta,T,5,0.5,(nchoosek(5,round(5/2))),2);
                        t_Gauss5(j,m) = toc;
                        % %                 %%%% Detect Julian
                        y =  squeeze(Yb(j,:,m));max_r =  Photons(j);
                        tic
                        [Z_Histo(j,m),W(j)] = detect(y(:),h,max_r);
                        t_Histo(j,m)=toc;
                    end
                    %                 pause
                end
                Z_Histo  = (Z_Histo>0);
                
                PFA_Gauss = mean((Z_Gauss==1).*(Z_GT==0),2);
                
                PFA_Gauss8 = mean((Z_Gauss8==1).*(Z_GT==0),2);
                
                PFA_Gauss5 = mean((Z_Gauss5==1).*(Z_GT==0),2);
                PFA_Histo = mean((Z_Histo==1).*(Z_GT==0),2);
                
                figure(2);
                subplot(2,1,Choice_Bshape)
                semilogx(Photons,mean(PFA_Histo,2),'r','LineWidth',1.5);hold on;
                semilogx(Photons,mean(PFA_Gauss,2),'bo-','LineWidth',1.5);hold on
                semilogx(Photons,mean(PFA_Gauss8,2),'g+-','LineWidth',1.5);
                %         semilogx(Photons,mean(PFA_Gauss5,2),'mo-','LineWidth',1.5);
                
                %         legend('ETD (M=10)','ETD (M=8)','ETD (M=5)','HTD')
                legend('HTD','ETD','GETD')
                ylabel('PFA '); xlabel('Photons');
                title(BackTitle{Choice_Bshape})
                set(gca,'FontSize',14);grid
                
%                 cd('./results')
%                 save(['Results_Synth_PFA_' BackTitle2{Choice_Bshape} '.mat'],'ImB','Yb','Z_GT','Photons','Z_Gauss','Z_Histo','t_Gauss','t_Histo','PFA_Gauss','PFA_Histo')
%                 saveName = ['ProbaFalseAlarm_HTD_ETD_GETD_Photons'];
%                 %                                 saveName = ['ProbaFalseAlarm_HTD_ETD_GETD_Photons_Sig6_T1000'];
%                 saveas(gcf, saveName, 'fig');
%                 print(saveName,'-dpng','-r512')
%                 cd ..
                
        end
    end
end

