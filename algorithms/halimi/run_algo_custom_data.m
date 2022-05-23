clearvars
useTargetDetect=0;
% F_gauss_sig_6=1,  F_real_proc=2
selirf=1;
run ../init.m
PPP = [0.1, 1, 5, 10, 50, 100];
SBR  = [0.1, 1, 5, 10, 50, 100];
T=K;

GuideD = 0;  %0 for Multiscale median, 1 for PC denoise
GuideI =  0; %0 for no filtering, 1 for DnCNN

dataDir = strcat('../../data_generation/',irfs(selirf));
outDir = strcat('../../results/halimi/',irfs(selirf));

if ~exist(outDir, 'dir')
    mkdir(outDir)
end

init_custom_data

%ppp=1; sbr=1;
LoopInd = 1;
n=1;
for k=1:2
    for ppp=PPP
        for sbr=SBR
            inFile=sprintf("Samples_%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f.mat", ...
                selectedScene, s_back{k}, K, downSam,ppp, sbr);
            path = strcat(dataDir,'/',inFile);
            if ~isfile(path)
                warning(strcat("Could not find file: ",inFile));
                %continue;
            end
            load(strcat(dataDir,'/',inFile),'Y');
            NbrePulse     = max(Y(:))*T ; %number of laser pulses, provided  in data or estimated as in here;
            Y=full(Y);
            Y=reshape(Y,row,col,[]);     

            fprintf('Processing case PPP=%i, SBR=%i', ppp, sbr);

            [Dscales, Rscales,   bT(LoopInd,:) ,  bN(:,LoopInd),sumB(LoopInd), Depth_CM] = Fct1_Extract_Estimates_withBack(Y,I_resol,F,IRFw,Tbin,Attack  ,trailing, mask);

            [zt,zt2,at22,zt3,Br,Eps_d, Eps_I] = Robust_Median_Bayesian(Dscales(:,:,:), (Rscales(:,:,:,:)),...
                I_resol,Neighbours_Med(:,1:indGraph_Med(2,2)),Neighbours_MedSR(:,1:indGraph_MedSR(2,2)) ,weights,local_shifts_Med, local_shifts_WM...
                ,Neighbours_WM, indGraph_WM ...
                ,ThreshDep,FramHist ,sigIRF,NbrePulse,Tbin,SR,GuideD,GuideI);       

            zt(isnan(zt))=0;
            Dest(:,:,LoopInd)    = reshape(zt,rowSR,colSR);
            Iest(:,:,LoopInd)  = reshape(zt2,[rowSR,colSR]);
            STD_d(:,:,LoopInd)   = reshape(Eps_d,rowSR,colSR);
            STD_I(:,:,LoopInd) = reshape(Eps_I,[rowSR,colSR]);
            dae=sum(sum(abs(Dest-Dref)))/(row*col);
            fprintf(', DAE=%f\n',dae);

            fprintf("Saving file %i/%i\n",n,length(PPP)*length(SBR)*2);
            outFile=sprintf("%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f.mat", selectedScene, s_back{k}, K, downSam, ppp, sbr); 
            save(strcat(outDir,'/',outFile), 'Dest', 'Iest', 'STD_d', 'STD_I', 'I_resol');   
            n=n+1;
        end 
    end
end
fprintf('Processing finished')