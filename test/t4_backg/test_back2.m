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

load ../../data_generation/Art_ref_img.mat
DOEs=zeros(length(PPP), length(SBR), 6);

k=2;
n=1;
for i=1:length(PPP)
    for j=1:length(SBR)
        ppp=PPP(i); sbr= SBR(j);

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
        fprintf("Processing case %i/%i\n",n,length(PPP)*length(SBR));

        [Dep, Refl] = estimateDepthTof(Y, neighboursSM, params, 0, 0);        
        DOEs(i,j,1) = mean(abs(Dep(:)-Dref_orig(:)));
        %[Dep, Refl] = estimateDepthTof(Y, neighboursSM, params, 1, 0);        
        %DOEs(i,j,2) = mean(abs(Dep(:)-Dref_orig(:)));
        for l=0:3
            [Dep, Refl] = estimateDepthHist(Y, F, neighboursSM, params, l);
            DOEs(i,j,l+3) = mean(abs(Dep(:)-Dref_orig(:)));
        end
        n=n+1;
    end
end

i=6;
ppp=PPP(i);
figure;plot(squeeze(DOEs(i,:,:)));
legend('TD wo. BackEst', 'TD with Backest', 'XCorr wo. BackEst', ...
    'XCrorr + BackEst Neighbours', 'XCorr + BackEst Profile','XCorr + BackEst Neighbours w. MXCorr');
xlabel("SBR"); xticklabels(SBR);
title(sprintf("DAE after cross correlation, PPP=%f", ppp));

DOEs_backup=DOEs;
save("back2_res","DOEs")