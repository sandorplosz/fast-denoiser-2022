function [Neighbours, indGraph,local_shifts0] =  Build_Graph_Neighbours_Array_v2(row,col,I_resol);
R        = length(I_resol);
N        = row*col;
Im       = reshape(1:N, row,col);
Fpix     = zeros(row,col,R);
SizeDup  = max(I_resol);
[ImExt]  = Duplicate_Borders(Im,SizeDup);
Neighbours = []; % Graph N+1  X  ndxSizeGraphs
indGraph = zeros(2,R); % Start/end for pixels
start     = 0;
local_shifts0 = [];
for r   =  1:R
    if(I_resol(r)==2)
        [ltemp ctemp]  = find([0 1 0;1 1 1;0 1 0]);
        local_shifts   = [ctemp ltemp]- 2;
    else
        [ltemp ctemp]  = find(ones(I_resol(r)));
        local_shifts   = [ctemp ltemp]- ceil(I_resol(r)/2);
    end
    
    if(r>1)
        local_shifts   = setdiff(local_shifts,local_shifts0,'rows');
    else
        nd = 1;
    end
    nd             = size(local_shifts,1);
     
    
    for i=1:nd
        ImExt2(:,:,i) = circshift(ImExt,[local_shifts(i,:)]);
    end
    
    
    %     Neighbours{r}  = reshape(ImExt2(SizeDup+1:end-SizeDup,SizeDup+1:end-SizeDup ,:),row*col,nd);
    Neighbours     = [Neighbours, ...
        [reshape(ImExt2(SizeDup+1:end-SizeDup,SizeDup+1:end-SizeDup ,:),row*col,nd); r*ones(1,nd)]];
    indGraph(:,r)  = [start+1; start+nd];start = start+nd;
    local_shifts0=[local_shifts0;local_shifts];
end