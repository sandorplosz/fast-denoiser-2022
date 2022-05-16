function [det, DiffLog, D_Xcorr,NbreGood] = detect_AH_NonUnifBack_v7_v2(tof,sigma2,Alpha,Beta,T,PhotStep,ProbPrior,limitC,varo, BackT, debug)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model  Unif + Gauss
%%% Limit the computations of c_k to the first similar elements
%%% truncate large combinations by considering indices closest to median
%% Input
%tof
% sigma2: width of IRF in bins^2
% Alpha,Beta: Prior parameters
% T:  length observation window
% PhotStep:  Nbre of photons to approximate computations (choose <12)
% ProbPrior: default 0.5
% limitC: related to PhotStep
% varo:   default 3
% BackT: a 1xT vector representing the background shape, if uniform Back=1/T*ones(1,T); ********
%% Output
% det: binary detection map, 
% DiffLog: difference log probas 
% D_Xcorr: cross-correlation depth 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


warning('off','all')
if length(tof)==0
    det = 0;DiffLog=-100; D_Xcorr=1; 
    Coeff=0;,aprev=0;  NbreGood =0;
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Study if TOF not empty
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Sort TOF: cluster first
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bord1     = min(3*sqrt(sigma2),round(T/4));%3*sqrt(sigma2);
    bord2     = max(T-3*sqrt(sigma2),round(3*T/4));%T-bord1;
    %logp1Fct = [ ];
    %logp0Fct = [ ];
    P         = length(tof); % Nbre photons
    %p1(1)     = ProbPrior;
    DevN2     = ceil(PhotStep/2);
    BackT     = log(BackT(tof)); %********
%     BackT     = BackT(:);        %********
    sumB      = sum(BackT);      %********
%     weight    = zeros(1,T);      %********
%     keyboard
%     weight([1:bord1, bord2:T])    = 30;   %********
%     weight    = BackT + weight;          %********
     
%     %%%%   matched filter (approximate!! )
     %v         = sort( exp(-  (tof - tof').^2/2/sigma2 - BackT(:)/2-BackT(:)'/2 ));       %********
     v= exp(- (tof - tof').^2/2/sigma2 - BackT(:)/2-BackT(:)'/2 );
     %[~, ind] = sort(sum(v),'descend');     %********
     [~, ind] = max(sum(v));     %********
     D_Xcorr   = tof(ind);     %********
      
     
%     %%%%  Log matched filter
%      v         = sort((tof - tof').^2);
%      [val ind] = sort(sum(v(1:min(P,PhotStep+1),:)));
%      D_Xcorr      = tof0(ind(1));%median(tof0);
    %%%%%
    
    [v2, ind]  = sort(abs(tof-D_Xcorr));
    tof      = tof(ind);
    BackT     = BackT(ind);     %********
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Compute Coefficients
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    B = betaln(Alpha,Beta);
    Xph    = 0:P-1;
    Coeff = betaln(Alpha+P-Xph,Beta+Xph) ... %- Xph*log(T) ...  from sum over i in page24 or p27b,28b     %********
        -0.5*log((P-Xph)) - ((P-Xph-1)/2)*log((2*pi*sigma2)) ... % from Cij in page 7     %********
        -B-log(T); % from constant in page 24 or p27b,28b     %********
%     Coeff(P+1) = betaln(Alpha,Beta+P)-log(T)*(P-1)-B-log(T);
    Coeff(P+1) = betaln(Alpha,Beta+P) -B + sumB;     %******** Is the only place where Back is used?
 
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Compute logp0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    logp0 = log(1-ProbPrior) + sumB;%- P*log(T);     %******** 
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Compute logp1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    is = min(P,round(1.2*PhotStep)):-1:1;
    vect=zeros(1,length(is));     
    for i=is
        %     for i=   0:P-1
        % %         nbreElements(i+1)  =nchoosek(P,P-i);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% Compute truncated Possib
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Logratio = 0;  
        %         if((nchoosek(P,P-i)) > 15000)% limitC )
        [Possib, cs]=getPossibMatrix2(P, i, PhotStep, limitC);
        NbreGood  = sum(v2 < ( max(v2(i),  varo*sqrt(sigma2) ) )) ;
        if cs==1            
            Logratio  =   max(log(nchoosek(max(i,NbreGood),i))-log(size(Possib,1)),0);
        end        
       
        if debug==1
        fprintf("i=%i. matrix size: %ix%i, NbreGood: %i, Logratio: %d\n", ...
            i,size(Possib,1),size(Possib,2), ...
            NbreGood, Logratio);        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% Compute log a:  aprev
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
        s      = tof(Possib);
        fb     = BackT(Possib);     %********
        if(size(s) ~= size(Possib))
            s = s';fb = fb';     %********
        end 
        fb2     = sumB - sum(fb,2);     %******** 
        TruncWind = mean(s,2);
        TruncWind = (TruncWind<bord1 | TruncWind>bord2 ); % to discard borders
        c = (Coeff(P-i+1))  - 30* TruncWind   + fb2 + (-0.5/sigma2* (sum(s.^2,2) - sum(s,2).^2 /i  ));     %********
          
        CsteMax     = max(c);
        %C{i+1}      = c;
        vect(i) = log(sum(exp(c -CsteMax))) + CsteMax  + Logratio;
        
    end
     
    aprev            = [vect,  Coeff(P+1)];
    CsteMax          = max(aprev);
    logp1            = log(ProbPrior) +  log(sum(exp(aprev -CsteMax))) + CsteMax;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%  Make decision
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if debug && NbreGood  ~= sum(v2 < ( varo*sqrt(sigma2) ) ) 
        fprintf("NbreGood not equal!\n");
    end
    DiffLog   =  ( logp1-logp0 );
    %det       = ( logp1-logp0 )> 0;
    det       = ( logp1-logp0 )>= 0;
        
end
