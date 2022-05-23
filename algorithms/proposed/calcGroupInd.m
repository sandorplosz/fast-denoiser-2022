gindexes=zeros(size(Tofmat,1),2);
numg=zeros(ceil(size(Tofmat,2)/32),1);
for i=1:size(Tofmat,1)
    t = Tofmat(i,:);
    t=t(t>1);
    if(~isempty(t))
        gi = floor((length(t)-1)/32);
        gindexes(i,1)= gi;
        gindexes(i,2)= numg(gi+1);        
        if i<20
            fprintf("Vector %i, size: %i, group: %i, index in group: %i\n", i-1, length(t), gi, numg(gi+1));
        end        
        numg(gi+1) = numg(gi+1) + 1;
    end
end