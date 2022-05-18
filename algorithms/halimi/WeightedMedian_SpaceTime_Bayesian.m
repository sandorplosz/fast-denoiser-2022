function [Dep vectMat weights0]  = WeightedMedian_SpaceTime_Bayesian(Dscale,Dtime,weights0,Eps_d,Neighbours_scale,Neighbours_time);

%%% Dscale scale images  r x c x L
%%% Dtime  past  images  r x c x T
%%% Eps_d  variance of size N x 1
%%% Wscale of size   N x L
%%% Wtime  of size  N x Ntime.T 
%%% Neighbours_scale of size N x Nscale
%%% Neighbours_time of size N x Ntime

%%% Note that temporal correlation is not used in this demo 


[N Nscale] = size(Neighbours_scale);
Ntime      = size(Neighbours_time,2);

[N,L]    = size(Dscale); 
Dmat1    = Dscale;% reshape(Dscale,N,L);
if(isempty(Dtime)==0)
    past           = size(Dtime,2);
    Dmat2          = Dtime;%reshape(Dtime,N,past); 
    weights0 = weights0(:,1:(Nscale *L+past*Ntime) );
else
    past = 0; 
    weights0 = weights0(:,1:Nscale *L);
end

%%% Weights
weights = weights0 ./sum(weights0,2);
Lw      = size(weights,2); 
weights = weights./ Eps_d(:); %N Lw
weights = weights';%Lw N

%  keyboard
%%% Images
vectMat        = [Dmat1(:,1)]; % N Lw 
vectMat        = [vectMat(Neighbours_scale)]; % N Lw
for t=2:L
    D2         = Dmat1(:,t);
    vectMat    = [vectMat, D2(Neighbours_scale)]; % N Lw
end 
for t=1:past
    D2         = Dmat2(:,t);
    vectMat    = [vectMat, D2(Neighbours_time)]; % N Lw
end 
 
%%% WMF
[vectOrderMat IndMat] = sort(vectMat,2); % N x T  
IndMat2               = IndMat' + repmat((0:Lw:(N-1)*Lw),Lw,1); % Lw N
WM2                   = weights(IndMat2(:)) ;
WM2                   = reshape(WM2,Lw,N )' .* (isnan(vectOrderMat)==0);%N Lw
WM2                   = WM2./sum(WM2,2); %N Lw
w0Mat                 = 1+sum(cumsum(WM2')<0.5); % 1xN
 
for n = 1:N
        Dep(n)  = vectOrderMat(n,w0Mat(n)) ;  % N x 1
end 