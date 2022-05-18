clearvars
useTargetDetect = 0; 
% F_gauss_sig_6=1,  F_real_proc=2
selirf=2;
init
estimateBackground=2;
PPP = [0.1, 1, 5, 10, 50, 100];
SBR  = [0.1, 1, 5, 10, 50, 100];
addpath('../../algorithms/proposed/')
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

ppp=10; sbr=10;
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
load ../../data_generation/Art_ref_img.mat

if useTargetDetect
    [Dep, Refl] = estimateDepthTof(Y, neighboursSM, params, estimateBackground);        
else
    [Dep, Refl] = estimateDepthHist(Y, F, neighboursSM, params, estimateBackground);
end

mean(abs(Dep(:)-Dref_orig(:))) % 66
a=Dep~=0;
mean(abs(Dep(a)-Dref_orig(a))) % 22.4
keyboard

% Load app results for comparison
%run mres_app.m

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

%sum(abs(Dest(:)-den_depth(:))>0.1)
%figure;imagesc(abs(Dest-den_depth))