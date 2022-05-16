function [F] = processF(F, T)
    F = F(:,100); F(F<(0.01*max(F(:))))=0;
    if T<length(F)
        h = [zeros(100,1); F(2:2:end); zeros(100,1)];
    else
        l = T-length(F);
        h = [zeros(floor(l/2),1); F; zeros(ceil(l/2),1)];
    end
    h  = h/sum(h);
    [~, MaxF] = max(h); 
    clear F;
    l = length(h);
    F=zeros(l,l);    
    for i=1:l
        F(:,i)  = circshift(h,[i-MaxF 0]);
    end
    F = F(1:T,1:T);
end