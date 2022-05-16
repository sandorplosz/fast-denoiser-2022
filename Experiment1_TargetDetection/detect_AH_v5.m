function [det, logp1, logp0,p1,p1Spat, logp1Fct, logp0Fct,aprev,C,MedT] = detect_AH_v5(tof,sigma2,Alpha,Beta,T,PhotStep,ProbPrior,limitC,varo);

%%% Model  Unif + Gauss
%%% Limit the computations of c_k to the first similar elements
%%% truncate large combinations by considering indices closest to median


warning('off','all')
if length(tof)==0
    det = 0;Prob1=0; Prob0=1;p1Spat=0;logp1=-100;logp0=0;logp1Fct = [ ];    p1(1)     = ProbPrior;aprev=0;C=0;MedT=0; 
    logp0Fct = [ ];
elseif(length(tof)==1)
    P         = 1;
    logp0     = log(1-ProbPrior) - P*log(T);
    logp1     = log(ProbPrior/T/beta(Alpha,Beta)   * ((beta(1+Alpha,Beta)+beta(Alpha,1+Beta))));
    bord1     = T/4;%3*sqrt(sigma2);
    bord2     = 3*T/4;%T-bord1;
%     if(tof<bord1 | tof>bord2) logp1=logp1-10; end
    if(abs(logp1 -logp0)>10^(-8) )
        det  = ( logp1-logp0 ) > 0;
    else
        det  = randi(2)-1;
    end
    Prob1=exp(logp1); Prob0=exp(logp0);p1Spat=Prob1;logp1Fct = [ ]; logp0Fct = [ ];aprev=0;C=0;MedT=tof;
    p1(1)     = ProbPrior;
     
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Study if TOF not empty
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Sort TOF: cluster first
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bord1     = T/4;%3*sqrt(sigma2);
    bord2     = 3*T/4;%T-bord1;
    logp1Fct = [ ];
    logp0Fct = [ ];
    P         = length(tof); % Nbre photons
    p1(1)     = ProbPrior;
    DevN2     = ceil(PhotStep/2);
    tof0      = tof;
    
%     %%%%   matched filter
%      v         = sort(1- exp(-(tof - tof').^2/2/sigma2));
%      [val ind] = sort(sum(v));
%      MedT      = tof0(ind(1));
%     %%%%  Log matched filter
     v         = sort((tof - tof').^2);
     [val ind] = sort(sum(v(1:min(P,PhotStep+1),:)));
     MedT      = tof0(ind(1));%median(tof0);
    %%%%%
      v         = sort( exp(-  (tof - tof').^2/2/sigma2   ));       %********
     [val ind] = sort(sum(v),'descend');     %********
     MedT   = tof0(ind(1));     %********
      
     
    
    [v2 ind]  = sort(abs(tof0-MedT));
    tof0      = tof0(ind); tof=tof0;
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Compute Coefficients
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    B = betaln(Alpha,Beta);
    Xph    = 0:P-1;
    Coeff = betaln(Alpha+P-Xph,Beta+Xph) - Xph*log(T) ...
        -0.5*log((P-Xph)) - ((P-Xph-1)/2)*log((2*pi*sigma2))-B-log(T);
    Coeff(P+1) = betaln(Alpha,Beta+P)-log(T)*(P-1)-B-log(T);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Compute logp0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    logp0 = log(1-ProbPrior) - P*log(T);
      
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Compute logp1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vect=[]; 
    for i= max(0,round(P-1.2*PhotStep)):P-1
        %     for i=   0:P-1
        % %         nbreElements(i+1)  =nchoosek(P,P-i);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% Compute truncated Possib
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Logratio = 0;  
        %         if((nchoosek(P,P-i)) > 15000)% limitC )
        if((nchoosek(P,P-i)) >  limitC ) %15000)%
             
 
            if((P-i)<= P/2  & (P-i)>DevN2)
                
                Possib = repmat((P-i):-1:1,limitC,1);
                for ii=1:min(limitC,P-i)
                    Possib(ii+1,ii) = P-i+1;
                end
                Possib = Possib(1:min(limitC,P-i),:);
                
            elseif((P-i)<= P/2  & (P-i)<=DevN2)
                Possib = nchoosek(1:min(22-2*(P-i),P),P-i);
                
            elseif((P-i)>P/2  & i>DevN2)
                vvv = nchoosek(1:PhotStep,DevN2);
                Possib =     [repmat(1:(P-i-DevN2),size(vvv,1),1)  (P-i-DevN2)+vvv ];
                
            elseif((P-i)>P/2  & i<=DevN2)
                k=PhotStep-i;vvv = nchoosek(1:PhotStep,PhotStep-i);
                Possib =     [repmat(1:(P-i-k),size(vvv,1),1)  (P-i-k)+vvv ];
            end
            
            
            %%% Correction coefficient for the approximation
            NbreGood  = sum(v2 <  ( max(v2(P-i),  varo*sqrt(sigma2) ) )) ;
            Logratio  =   max(log(nchoosek(max(P-i,NbreGood),P-i))-log(size(Possib,1)),0);
            
            
        else
            Possib = nchoosek(1:P,P-i); sss(i+1,1:2)=size(Possib);
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% Compute log a:  aprev
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
        s      = tof(Possib);
        if(size(s) ~= size(Possib))
            s = s';
        end
        TruncWind = mean(s,2);
        TruncWind = (TruncWind<bord1 | TruncWind>bord2 ); % to discard borders
        c = (Coeff(i+1))  - 30* TruncWind  + (-0.5/sigma2* (sum(s.^2,2) - sum(s,2).^2 /(P-i)  ));
        CsteMax     = max(c);
        C{i+1}      = c;
%         aprev(i+1)  = log(sum(exp(c -CsteMax))) + CsteMax  + Logratio; % in log scale
        vect = [vect, log(sum(exp(c -CsteMax))) + CsteMax  + Logratio];
        
    end
%     aprev            = setdiff(aprev,0);
%     aprev(end+1)     = Coeff(P+1);  % 1 x P+1  in log scale
    aprev            = [vect,  Coeff(P+1)];
    CsteMax          = max(aprev);
    logp1            = log(ProbPrior) +  log(sum(exp(aprev -CsteMax))) + CsteMax;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%  Make decision
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    c         = max(min(100,exp(logp1-logp0)), 0.01); %  0.01< p1/p0 <100
    p1Spat    =  c/(1+c);%
    logp1Fct  = [logp1Fct, logp1];
    logp0Fct  = [logp0Fct, logp0];
    if(logp1 ~= logp0)
        det  = ( logp1-logp0 ) > 0; 
    else
        det  = randi(2)-1;
    end
end








 