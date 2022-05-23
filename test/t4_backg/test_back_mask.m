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

ppp=100; sbr=1;
k=2;

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

%%%%%%
disp('    Creation TOF from histogram')
n = 1;
N=row*col;
Nmax = max(max(sum(Y,3)));
for j=   1:col
    for i=   1:row
        y = squeeze(Y(i,j,:));
        Tof=[];
        ind = find(y>0);
        for l=1:length(ind)
            Tof = [Tof kron(ind(l),ones(1,y(ind(l))) )];
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

DD=zeros(1,N);
Refl=zeros(1,N);
Dep=zeros(1,N);
Start=tic;
for n= 1:N
    TOFn  = Tofmat(n,:);
    TOFn  = TOFn(TOFn>1);
    [~, ~, DD(n), Refl(n)] = ...
        detect_AH_NonUnifBack_v7_v2(TOFn, params.sigIRF^2, params.Alpha,params.Beta,K, ...
        params.NbrePhoton, params.ProbPrior(n), params.limitC, ...
        3, 1/K*ones(1,K), 0);%5)
end
toc(Start);

%dd=reshape(DD,row,col);
%figure;imagesc(dd);
fprintf("DAE after XCorr: %f\n", log10(mean(abs(DD(:)-Dref_orig(:)))*params.Tbin *3*10^8/2));

Back = 1/K*ones(1,K);
DD      = DD(:);
Dist    = abs(DD - DD(Neighbours.neighb(:,1:Neighbours.indGraph(2,2))));
Dist2   = sort(Dist,2,'ascend');
mask    = sum(Dist2(:,2:4),2)>150;   %Dist2(:,2)>5;
fprintf("Sum mask: %i\n", sum(mask));
%if(sum(mask)>min(2*K,0.1*N))
    Back    = hist(DD.*mask(:),1:K);
    Back(1) = Back(2);
    Back    = movmean(Back,[params.Attack, params.trailing]);               
    Back    = Back/sum(Back);
    figure;plot(Back);

    Start=tic;
    for n= 1:N
        TOFn  = Tofmat(n,:);
        TOFn  = TOFn(TOFn>1);
        [~, ~, Dep(n), Refl(n)] = ...
            detect_AH_NonUnifBack_v7_v2(TOFn, params.sigIRF^2, params.Alpha,params.Beta,K, ...
            params.NbrePhoton, params.ProbPrior(n), params.limitC, ...
            3, Back, 0);%5)
    end
    fprintf("DAE after Back est: %f\n", log10(mean(abs(Dep(:)-Dref_orig(:)))*params.Tbin *3*10^8/2));
    toc(Start);    
% else
%     Dep=DD(:);
%     fprintf("Background kept uniform!\n");
% end

% Uniform Backg
%--------------
% PPP=50, SBR=0.1
% DAE after XCorr= -0.376, 
% DAE after den.=  -1.098  
% Threshold: 80, sumMask=42179, after Back est: -0.13, after den: -0.39   
% Threshold: 150, sumMask=37134, after Back est: -0.48, after den: -1.4   
% Threshold: 200, sumMask=33624, after Back est: -0.85, after den: -1.57

% PPP=10, SBR=0.1
% DAE after XCorr= -0.076581
% DAE after den.=  -0.307380
% Threshold: 80, sumMask=80366, after Back est: -0.141903, after den: -0.484707   
% Threshold: 200, sumMask= 63544, after Back est: -0.296503, after den: -0.993630
% Threshold: 250, sumMask= 55917, after Back est: -0.372948, after den: -1.080028
% Threshold: 300, sumMask= 48137, after Back est: -0.451927, after den: -1.077641

% PPP=0.1, SBR=0.1
% DAE after XCorr= 0.348937
% DAE after den.=  -0.089748
% Threshold: 50, sumMask=9047 , after Back est: 0.3489, after den: -0.102
% Threshold: 100, sumMask=8909 , after Back est: 0.3487, after den: -0.102
% Threshold: 150, sumMask=8766, after Back est: 0.348713, after den: -0.102054
% Threshold: 200, sumMask=8606 , after Back est: same, after den: same

% Gamma backg
%------------
% PPP=0.1, SBR=0.1
% DAE after XCorr= 0.36
% DAE after den.= 0.21
% Threshold: 50, sumMask=8873 , after Back est: 0.36, after den: 0.20
% Threshold: 100, sumMask=8664 , after Back est: 0.36, after den: 0.20
% Threshold: 200, sumMask=7913 , after Back est: same, after den: same
% Threshold: 300, sumMask=6946 , after Back est: same, after den: same

% PPP=10, SBR=0.1
% DAE after XCorr= 0.205
% DAE after den.= 0.23
% Threshold: 50, sumMask=83911 , after Back est: -0.11, after den: -0.11 
% Threshold: 100, sumMask=60339 , after Back est: -0.09 , after den: -0.1 
% Threshold: 150, sumMask=41567 , after Back est: -0.02, after den: -0.06 
% Threshold: 200, sumMask=30122 , after Back est: 0.1, after den: 0.12
% Threshold: 300, sumMask=18797 , after Back est: 0.22, after den: 0.28

% PPP=100, SBR=10
% DAE after XCorr= -2.31
% DAE after den.= -2.34
% Threshold: 150, sumMask= 18, after Back est: -1.99, after den: -2.13

% PPP=100, SBR=1
% DAE after XCorr= -1.7
% DAE after den.= -2.2
% Threshold: 100, sumMask= 683, after Back est: -0.35, after den: -0.44
% Threshold: 150, sumMask=606, after Back est: -0.67, after den: -0.92 

% PPP=100, SBR=0.1
% DAE after XCorr= 0.076
% DAE after den.= 0.024
% Threshold: 100, sumMask=31167 , after Back est: -0.17 , after den: -0.13
% Threshold: 150, sumMask=21147 , after Back est: -0.17, after den: -0.13

% Threshold: 00, sumMask= , after Back est: , after den: 

%%%%%%

% if useTargetDetect
%     [Dep, Refl] = estimateDepthTof(Y, neighboursSM, params, estimateBackground, 0);        
% else
%     [Dep, Refl] = estimateDepthHist(Y, F, neighboursSM, params, estimateBackground);
% end

if(0)   
    % 197.9
    a=Dep~=0;
    mean(abs(Dep(a)-Dref_orig(a))) % 22.4
    keyboard
    
    % Load app results for comparison
    run mres_app.m
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

fprintf("DAE after denoising: %f", log10(mean(abs(Dest(:)-Dref(:)))));
%log10(mean(abs(Dest(:)-Dref(:))))
%sum(abs(Dest(:)-den_depth(:))>0.1)
%figure;imagesc(abs(Dest-den_depth))