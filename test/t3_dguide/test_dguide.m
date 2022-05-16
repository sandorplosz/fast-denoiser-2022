clearvars
Background = 0;  %0 for uniform backgroun, 1 Gamma-shaped background
useTargetDetect = 0;
estimateBackground = 0;
init

outDir = '/home/ps2014/Development/matlab/fast_denoiser_2022/results_algo';

addpath('/home/ps2014/Development/matlab/fast_denoiser_2022/algo')
load /home/ps2014/Development/matlab/fast_denoiser_2022/data_generation/F_real2_100s
F = processF(F,K);

Neighbours.I_resol = [1 5 7 9] ;  % size of spatial correlations Requires
%Neighbours.I_resol = [1 7 13] ;  % size of spatial correlations Requires
[Neighbours.neighb, Neighbours.indGraph, Neighbours.local_shifts]  =  Build_Graph_Neighbours_Array_v2(row,col,[1 Neighbours.I_resol(2:end)]); % Define graph of correlations between pixels
Neighbours.neighb  = Neighbours.neighb(1:end-1,:);
neighboursSM = Neighbours.neighb(:,1:Neighbours.indGraph(2,3));
s_filters="";
for i=1:length(Neighbours.I_resol)-1
    s_filters=strcat(s_filters, int2str(Neighbours.I_resol(i)), ", ");
end
s_filters=strcat(s_filters, int2str(Neighbours.I_resol(end)));

%ppp=10;
%sbr=10;

% DGuide selection
% 1. Select first
% 2. Select median
% 3. Select mean
% 4. Select mean2:end
DAE = zeros(length(PPP),length(SBR), 4);
IAE = zeros(length(PPP),length(SBR), 4);

n=1;
for i=1:length(PPP)
    for j=1:length(SBR)
        ppp=PPP(i); sbr= SBR(j);
        inFile=sprintf("%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f.mat", ...
            selectedScene, backName, K, downSam,ppp, sbr);
        path = strcat(dataDir,'/',inFile);
        if ~isfile(path)
            warning(strcat("Could not find file: ",inFile));
            %continue;
        end
        load(strcat(dataDir,'/',inFile),'Y');
        Y=full(Y);
        Y=reshape(Y,row,col,[]);
        fprintf("Processing file %i/%i\n",n,length(PPP)*length(SBR));
        %VectY = squeeze(sum(sum(Y,1),2)); sB = mean(VectY(1:300))*1024; sY=sum(VectY); sbrR=(sY-sB)/sB, pppR = sY/row/col
        %figure; plot(VectY);
        
        if useTargetDetect
            [Dep, Refl] = estimateDepthTof(Y, neighboursSM, params, estimateBackground);        
        else
            [Dep, Refl] = estimateDepthHist(Y, F, neighboursSM, params, estimateBackground);
        end
        Dscales = Dep * params.Tbin *3*10^8/2;
        
        fprintf('Running denoise.');        
        Dtype = 'Estimates';%'Hist_Back';        
        Start=tic;
        for k=1:4
            [Dest, Rest, ~, D_uncert, R_uncert] = Fct_Denoise_WMF_v4(Dep,Refl,Neighbours.I_resol,F,params.IRFw, params.Tbin, ...
                params.Attack, params.trailing, Neighbours.neighb, Neighbours.indGraph, Neighbours.local_shifts,Dtype, params.ThreshDep,0,0,params.CvPC, k);
            dae=sum(sum(abs(Dest-Dref)))/(row*col);
            DAE(i,j,k)=dae;
            IAE(i,j,k)=sum(sum(abs(Rest-Iref)))/(row*col);
            fprintf('.');            
        end
        toc(Start);        
        n=n+1;
    end 
end

s1={'select 1st', 'select median', 'select mean', 'select mean2'};
ires = Neighbours.I_resol;
outfname = sprintf("error_%s_%s_%s",selectedScene,backName,s_backest);
save(strcat(outfname,'.m'),'ires', 'DAE', 'IAE', 'PPP', 'SBR', 's1')

figure; sgtitle(sprintf("%s, %s, filters=[%s], estimateBackground=%i",selectedScene, s_back,s_filters,estimateBackground));
for k=1:4
    subplot(1,4,k);
    contourf(log10(PPP),log10(SBR),log10(DAE(:,:,k)));
    ylabel('SBR');xlabel('Photons');caxis([-3 0.5])
    title(strcat("Dguide ",s1{k}));
    xticks([-1 0 1 2]);xticklabels({'10^{-1}','1','10^1','10^2'})
    yticks([-1 0 1 2]);yticklabels({'10^{-1}','1','10^1','10^2'})
    colorbar
end
saveas(gcf,strcat(outfname,'.fig'));
saveas(gcf,strcat(outfname,'.png'));