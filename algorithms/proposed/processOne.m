clearvars
useTargetDetect = 1; 
% F_gauss_sig_6=1,  F_real_proc=2
selirf=1;
init
estimateBackground=1;
PPP = [0.1, 1, 5, 10, 50, 100];
SBR  = [0.1, 1, 5, 10, 50, 100];
addpath('../../algorithms/proposed/')
addpath ~/Development/Algorithm2_Parallel_TargetDetection_AH/build_cuda/
%addpath('/home/ps2014/Development/Algorithm2_Parallel_TargetDetection_AH/build_cuda')

Neighbours.I_resol = [1 3 7] ;  % size of spatial correlations Requires
[Neighbours.neighb, Neighbours.indGraph, Neighbours.local_shifts]  =  Build_Graph_Neighbours_Array_v2(row,col,[1 Neighbours.I_resol(2:end)]); % Define graph of correlations between pixels
Neighbours.neighb  = Neighbours.neighb(1:end-1,:);
neighboursSM = Neighbours.neighb(:,1:Neighbours.indGraph(2,2));
s_filters="";
for i=1:length(Neighbours.I_resol)-1
    s_filters=strcat(s_filters, int2str(Neighbours.I_resol(i)), ", ");
end
s_filters=strcat(s_filters, int2str(Neighbours.I_resol(end)));

ppp=1; sbr=10;
k=1;
Lev_S = sbr*ppp/(1+sbr);
Iref = IrefGray * Lev_S;

inFile=sprintf("Samples_%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f.mat", ...
        selectedScene, s_back{k}, K, downSam,ppp, sbr);
path = strcat(dataDir,'/',inFile);
if ~isfile(path)
    warning(strcat("Could not find file: ",inFile));
    %continue;
end
load(strcat(dataDir,'/',inFile),'Y');
Y=full(Y);
Y=reshape(Y,row,col,[]);

[row, col, T] = size(Y);
N=row*col;

disp('    Creation TOF from histogram')
n = 1;
Nmax = max(max(sum(Y,3)));
for j=   1:col
    for i=   1:row
        y = squeeze(Y(i,j,:));
        Tof=[];
        ind = find(y>0);
        for k=1:length(ind)
            Tof = [Tof kron(ind(k),ones(1,y(ind(k))) )];
        end
        Tof    = Tof(randperm(length(Tof)));            
        tof{n} = Tof;
        n      = n+1;
    end
end
Tofmat = zeros(N,Nmax);
for n = 1:N
    %ind = randperm(Nmax,nt(n));
    Tofmat(n,1:length(tof{n})) = sort(tof{n});
end

p=params;
DD=zeros(1,N);rr=zeros(1,N);
Start=tic;
for n= 1:N
    TOFn  = Tofmat(n,:);
    TOFn  = TOFn(TOFn>1);
    [~, ~, DD(n), rr(n)] = ...
        detect_AH_NonUnifBack_v7_v2(TOFn   ,p.sigIRF^2,p.Alpha,p.Beta,T,p.NbrePhoton,p.ProbPrior(n),p.limitC, 3, 1/T*ones(1,T), 0);%5)
end

DD      = DD(:);
Dist    = abs(DD - DD(neighboursSM));
Dist2   = sort(Dist,2,'ascend');
mask    = sum(Dist2(:,2:4),2)>150;   %Dist2(:,2)>5;
if(sum(mask)>min(2*T,0.1*N))
    Back    = hist(DD.*mask(:),1:T);
    Back(1) = Back(2);
    %Back    = movmean(Back,[1,2]);
    Back    = movmean(Back,[p.Attack, p.trailing]);               
    Back    = Back/sum(Back);
end

if(1)
    run mres_td.m
    td_dep1(td_dep1==0)=1;
    sum(DD~=td_dep1(:)') 
    a=find(DD(:)~=td_dep1(:));
    i=1; n=a(i); [n, DD(n), td_dep1(n)]
    
    figure;plot(abs(Back(:)-td_backg(:)))
    a=[Back(:), td_backg(:)];
    Back=td_backg(:);
end

Det_TD=zeros(1,N); DiffLogP=zeros(1,N);
Dep=zeros(1,N); Refl=zeros(1,N);
disp('    Running detection')
Start=tic;
for n= 1:N                    
    TOFn  = Tofmat(n,:);%(Iter-2)*Group+1:(Iter-1)*Group);
    TOFn  = TOFn(TOFn>1);
    TOFn=sort(TOFn,'descend');
    [Det_TD(n), DiffLogP(n), Dep(n),Refl(n)] = ...
        detect_AH_NonUnifBack_v7_v2(TOFn   ,p.sigIRF^2,p.Alpha,p.Beta,T,p.NbrePhoton,p.ProbPrior(n),p.limitC, 3, Back,0);%5)
end

if(1)  
    sum(Refl(:)~=td_refl(:))  
    b=find(rr~=td_refl(:)');    
    i=1; n=b(i); [n, rr(n), td_refl(n)]
    TOFn  = Tofmat(n,:);
    TOFn  = TOFn(TOFn>1);
    %DD=td_dep1(:)';
    fprintf("DAE_app after XCorr: %f\n", log10(mean(abs(td_dep2(:)-Dref_orig(:)))*params.Tbin *3*10^8/2));
    fprintf("IAE_app after XCorr: %f\n", log10(mean(abs(td_refl(:)-Iref(:)))/mean(Iref(:))));
    td_difflog(td_difflog==0)=-100;
    c=td_difflog(:)>0; d=DiffLogP(:)>0;sum(c~=d)
end

dopostproc=0;
if(dopostproc)
        indd=(sum(Y(:,:,2:end),3) == 1);  % AH-mod
        Det_TD(indd) = 1;            % AH-mod    
        Dep  = Dep(:).*Det_TD(:);
        ind     = (Det_TD(:)==0);
        r0      = reshape(Y(:,:,1),N,1);
        %R_TD(ind)    = r0(ind);
        Refl(ind)    = 10-r0(ind); Refl(r0>=10)=1; % AH-mod    
end


% if useTargetDetect
%     [Dep, Refl] = estimateDepthTof(Y, neighboursSM, params, estimateBackground, 0);        
% else
%     [Dep, Refl] = estimateDepthHist(Y, F, neighboursSM, params, estimateBackground);
% end

fprintf("DAE after XCorr: %f\n", log10(mean(abs(Dep(:)-Dref_orig(:)))*params.Tbin *3*10^8/2));
fprintf("IAE after XCorr: %f\n", log10(mean(abs(Refl(:)-Iref(:)))/mean(Iref(:))));

if(0)   
    % Load app results for comparison
    run mres_app.m
    sum(Dep(:)~=td_depth(:))
    fprintf("Algo DAE after XCorr: %f\n", log10(mean(abs(td_depth(:)-Dref_orig(:)))*params.Tbin *3*10^8/2));

    % 197.9
    a=Dep~=0;
    mean(abs(Dep(a)-Dref_orig(a))) % 22.4
    keyboard
    

    %run ../../results/proposed/F_real_proc/Art_UnifBack_K_1024_DownS_2_PPP_50_000_SBR_10_000_app.m
    log10(mean(abs(td_depth(:)-Dref_orig(:)))*params.Tbin *3*10^8/2)
    log10(mean(abs(den_depth(:)-Dref(:))))
end

%Dep=td_depth;Refl=td_refl;

if(0)
    % Check what causes the differences
    Yt=reshape(Y,row*col,[]);
    sum(Dep(:)~=td_depth(:))
    da=find(Dep~=td_depth);
    db=(abs(Dep(da)-td_depth(da)));
    for i=1:length(da)
        c=da(i);
        %[Dep(c),td_depth(c)]    
        y=find(Yt(c,:)); 
        % Mosly we see 2 photons, and different ones are selected, this is normal
        % In other cases the crossCorr maximum val is the same for the
        % other selected value, see:
        % [mv,mi]=sort(F*Yt(c,:)','descend');
        if(length(y)>2)
            fprintf("%i. diff!\n",i);
            break
        end
    end    
    sum(Refl(:)~=td_refl(:))
    ra=find(Refl~=td_refl);
    rb=abs(Refl(ra)-td_refl(ra));
    c=ra(7);
    [Refl(c), td_refl(c)]
    [Dep(c),td_depth(c)]
end

% Continue with the app output
%Dep=td_depth; Refl=td_refl;

Dtype = 'Estimates';
[Dest, Rest, ~, D_uncert, R_uncert] = Fct_Denoise_WMF_v3(Dep,Refl,Neighbours, params,F,Dtype, 0,0,params.CvPC, 3, 1);

fprintf("DAE after denoising: %f\n", log10(mean(abs(Dest(:)-Dref(:)))));
fprintf("IAE after denoising: %f\n", log10(mean(abs(Rest(:)-Iref(:)))/mean(Iref(:))));
%log10(mean(abs(Dest(:)-Dref(:))))
%sum(abs(Dest(:)-den_depth(:))>0.1)
%figure;imagesc(abs(Dest-den_depth))