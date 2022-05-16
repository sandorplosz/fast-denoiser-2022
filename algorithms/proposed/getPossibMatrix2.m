function [Possib, cs] = getPossibMatrix2(P, ind, PhotStep, limitC)

DevN2 = ceil(PhotStep/2);
Possib = 0;
s=round(1.2*PhotStep);

%is = min(P,s):-1:1;
if ind > s
    return;
end
    
%i=is(ind);
i=ind;
if((nchoosek(P,i)) >  limitC ) %15000)   % nchoosek(76,12...1)           
    if((i)<= P/2  && (i)>DevN2)              
        %fprintf("Case1. P: %i, i: %i\n",P,i);
        Possib = repmat((i):-1:1,limitC,1); % Matrix sor:(12...1) x 252 ismetelve
        for ii=1:min(limitC,i) % 1:12
            Possib(ii+1,ii) = i+1; % 2..13.sor 1..12. eleme = 13..2
        end
        Possib = Possib(1:min(limitC,i),:); % 1:12.sor

    elseif(i<= P/2  && i<=DevN2)        
        %fprintf("Case2. P: %i, i: %i\n",P,i);
        Possib = nchoosek(1:min(22-2*(i),P),i);
    elseif(i>P/2  && (P-i)>DevN2)        
        %fprintf("Case3. P: %i, i: %i\n",P,i);
        vvv = nchoosek(1:PhotStep,DevN2);
        Possib =     [repmat(1:(i-DevN2),size(vvv,1),1)  (i-DevN2)+vvv ];

    elseif(i>P/2  && (P-i)<=DevN2)                
        k=PhotStep-P+i;vvv = nchoosek(1:PhotStep,k);
        %fprintf("Case4. P: %i, i: %i, k: %i\n",P,i,k);
        Possib =     [repmat(1:(i-k),size(vvv,1),1)  (i-k)+vvv ];
    end               
    cs=1;
else
    %fprintf("Case5. P: %i, i: %i\n",P,i);
    Possib = nchoosek(1:P,i);
    cs=0;
end        
end