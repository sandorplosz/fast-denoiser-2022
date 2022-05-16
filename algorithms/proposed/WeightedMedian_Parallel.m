function Dep   = WeightedMedian_Parallel(D, D2, N,weights,Neighbours_WM, indGraph_WM,Filt);
%%% D current image r x c
%%% D2 past images  r x c x t
%%% Window size of window
%%% weightCurr: default 0.3
%%% N number of pixels
[row,col]      = size(D);
past           = 0;
Lw             = size(weights,2);
% weightsMat     = repmat(weights,N,1)';
weightsMat     = weights'; 
  
D              = D(:);
Dmat2          = reshape(D2,N,past);
vectMat        = [D(Neighbours_WM)]; % N t
for t=1:past
    D2            =  Dmat2(:,t);
    vectMat       = [vectMat, D2(Neighbours_WM)]; % N t
end
 
[vectOrderMat IndMat] = sort(vectMat,2); % N x 1  
IndMat2               = IndMat' + repmat((0:Lw:(N-1)*Lw),Lw,1);
WM2                   = weightsMat(IndMat2(:)) ;
WM2                   = reshape(WM2,Lw,row*col )';
if(Filt==1) % Weighted median
    w0Mat                 = 1+sum(cumsum(WM2')<0.5); % 1xN
    n=1;
    for c= 1:col
        for r= 1:row
            Dep(r,c)  = vectOrderMat(n,w0Mat(n)) ;  % N x 1
            n=n+1;
        end
    end
elseif(Filt==0) % Weighted averaging
    Dep = sum( vectOrderMat .* WM2  ,2);
    Dep = reshape(Dep,row,col);
     
end






% % function Dep   = WeightedMedian_SpaceTime(D, D2, N,weights,local_shifts,Neighbours_WM, indGraph_WM);
% % %%% D current image r x c
% % %%% D2 past images  r x c x t
% % %%% Window size of window
% % %%% weightCurr: default 0.3
% % %%% N number of pixels
% % [row,col]      = size(D);
% % past           = size(D2,3);
% % % Window=5;
% % % [ltemp ctemp]  = find(ones((Window)));
% % % local_shifts   = [ctemp ltemp]- ceil(Window/2);
% % nd             = size(local_shifts,1);
% % 
% % % nd             = size(local_shifts,1);
% % % p              = weightCurr;
% % % weights        = [p (1-p)/(2*nd-1)*ones(1,(2*nd-1))];
% % % weights        = [1./sum(abs(local_shifts(setdiff(1:nd,ceil(nd/2)),:).^2)') ];
% % % for t=1:past
% % %     weights = [weights 1./(t+sum(abs(local_shifts.^2)'))];
% % % end
% % % weights    = [p weights/sum(weights(:))*(1-p)];
% % Lw         = length(weights);
% % weightsMat = repmat(weights,N,1)';
% %   
% % D             = D(:);
% % Dmat2         = reshape(D2,N,past);
% % vectMat       = [D(Neighbours_WM)]; % N t
% % for t=1:past
% %     D2            =  Dmat2(:,t);
% %     vectMat       = [vectMat, D2(Neighbours_WM)]; % N t
% % end
% % % dtemp(:,:,1)  = D; 
% %  
% % 
% % % pp = 2;
% % % for v=[1:floor(nd/2) floor(nd/2+2):nd]
% % %     dtemp(:,:,pp) =circshift(dtemp(:,:,1),[local_shifts(v,:)]);
% % %     pp=pp+1;
% % % end;
% % % for t=1:past
% % %     for v=1:nd
% % %         dtemp(:,:,pp) =circshift(D2(:,:,1,t),[local_shifts(v,:)]);
% % %         pp=pp+1;
% % %     end;
% % % end
% % % vectMat               = reshape(dtemp,[size(dtemp,1)*size(dtemp,2),size(dtemp,3)]);% N t
% % 
% % [vectOrderMat IndMat] = sort(vectMat,2); % N x 1
% % % IndMat2               = IndMat+repmat((0:Lw:(N-1)*Lw)',1,Lw); % N 50
% % % IndMat2               = IndMat2';
% % IndMat2               = IndMat' + repmat((0:Lw:(N-1)*Lw),Lw,1);
% % WM2                   = weightsMat(IndMat2(:)) ;
% % WM2                   = reshape(WM2,Lw,row*col )';
% % w0Mat                 = 1+sum(cumsum(WM2')<0.5); % 1xN
% % n=1;
% % for c= 1:col
% %     for r= 1:row
% %         Dep(r,c)  = vectOrderMat(n,w0Mat(n)) ;  % N x 1
% %         n=n+1;
% %     end
% % end
% %  
% % 
% %  