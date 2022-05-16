clearvars
addpath('/home/ps2014/Development/matlab/fast_denoiser_2022/algo')
Background = 1;  %0 for uniform backgroun, 1 Gamma-shaped background
useTargetDetect = 0;
init

PPP = [0.1, 0.5, 1, 5, 10];
SBR  = [0.1, 0.5, 1, 5, 10];

outDir = '/home/ps2014/Development/matlab/fast_denoiser_2022/results_algo';
load /home/ps2014/Development/matlab/fast_denoiser_2022/data_generation/F_real2_100s
F = processF(F,K);

Neighbours.I_resol = [1 3 5] ;  % size of spatial correlations Requires
%Neighbours.I_resol = [1 7 13] ;  % size of spatial correlations Requires
[Neighbours.neighb, Neighbours.indGraph, Neighbours.local_shifts]  =  Build_Graph_Neighbours_Array_v2(row,col,[1 Neighbours.I_resol(2:end)]); % Define graph of correlations between pixels
Neighbours.neighb  = Neighbours.neighb(1:end-1,:);
neighboursSM = Neighbours.neighb(:,1:Neighbours.indGraph(2,2));
s_filters=""; s_filters2="";
for i=1:length(Neighbours.I_resol)-1
    s_filters=strcat(s_filters, int2str(Neighbours.I_resol(i)), ", ");
    s_filters2=strcat(s_filters2, int2str(Neighbours.I_resol(i)), "_");
end
s_filters=strcat(s_filters, int2str(Neighbours.I_resol(end)));
s_filters2=strcat(s_filters2, int2str(Neighbours.I_resol(end)));

%ppp=10;sbr=0.1;

DAE = zeros(length(PPP),length(SBR), 2, 3);
DAE_TD = zeros(length(PPP),length(SBR), 2, 3);
IAE = zeros(length(PPP),length(SBR), 2, 3);

n=1;    
for i=1:length(PPP)
    for j=1:length(SBR)
        for k=1:2 %1 for uniform backgroun, 2 Gamma-shaped background
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
            %Ysum=squeeze(sum(sum(Y,1),2)); figure;plot(Ysum);
            fprintf("Processing PPP=%f, SBR=%f (%i/%i)\n",PPP(i), SBR(j), n,numel(DAE));

            for l=1:3 % estimateBackground
                if useTargetDetect
                    [Dep, Refl] = estimateDepthTof(Y, neighboursSM, params, l);        
                else
                    [Dep, Refl] = estimateDepthHist(Y, F, Neighbours.neighb, params, l, 0, Dref);
                end
                Dscales = Dep * params.Tbin *3*10^8/2;
                dae=sum(sum(abs(Dscales-Dref)))/(row*col);                
                fprintf("DAE after est: %f\n", dae);
                DAE_TD(i,j,k,l)=dae;
                fprintf('Running denoise...\n');        
                Dtype = 'Estimates';%'Hist_Back';        
                Start=tic;
                [Dest, Rest, ~, D_uncert, R_uncert] = Fct_Denoise_WMF_v3(Dep,Refl,Neighbours,params,F,Dtype, 0,0,params.CvPC, 3, 1);
                dae=sum(sum(abs(Dest-Dref)))/(row*col);
                fprintf("DAE after denoise: %f\n", dae);
                DAE(i,j,k,l)=dae;
                IAE(i,j,k,l)=sum(sum(abs(Rest-IrefGray)))/(row*col);
                fprintf('.');
                toc(Start);        
                n=n+1;
            end 
        end
    end
end

s1 = {'Uniform backg.', 'Gamma backg.'};
s2={'No backg. estimation', 'Backg. estimation 1', 'Backg. estimation 2', 'Backg. estimation 3'};
if 0

end

figure; sgtitle(sprintf("Test with different background and estimation methods, %s, filters=[%s]",selectedScene, s_filters));
n=1;
for k=1:size(DAE,3)
    for l=1:size(DAE,4)
        subplot(size(DAE,3),size(DAE,4),n);
        contourf(log10(PPP),log10(SBR),log10(squeeze(DAE(:,:,k,l)')));
        ylabel('SBR');xlabel('Photons');caxis([-3 0.5])
        title(strcat(s1{k}, " ", s2{l}));
        %xticks(1:length(PPP));xticklabels(string(PPP))
        %yticks(1:length(SBR));yticklabels(string(SBR))
        xticks([-1 0 1]);xticklabels({'10^{-1}','10^0','10^1'})
        yticks([-1 0 1 2]);yticklabels({'10^{-1}','1','10^1','10^2'})
        colorbar
        n=n+1;
    end
end
if(1)
    ires = Neighbours.I_resol;
    outfname = sprintf("error_%s_backest_%s",selectedScene,s_filters2);
    save(strcat(outfname,'.m'),'ires', 'DAE', 'IAE', 'PPP', 'SBR', 's1', 's2')    
    saveas(gcf,strcat(outfname,'.fig'));
end