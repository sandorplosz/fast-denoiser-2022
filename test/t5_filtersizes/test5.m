clearvars
useTargetDetect = 0; 
init
estimateBackground=2;
PPP = [0.1, 0.5, 1, 5, 10];
SBR  = [0.1, 0.5, 1, 5, 10];

outDir = '/home/ps2014/Development/matlab/fast_denoiser_2022/results_algo';
addpath('/home/ps2014/Development/matlab/fast_denoiser_2022/algo')

filtersizes = { ...
    [1 3 5], [1 3 7], [1 3 9], [1 5 9], [1 7 11], [1 9 13], [1 5 9 11], [1 9 13 15]};

DAE = zeros(length(PPP),length(SBR),2);
IAE = zeros(length(PPP),length(SBR),2);

n=1;
for fs=filtersizes
    Neighbours.I_resol = fs{:};
    [Neighbours.neighb, Neighbours.indGraph, Neighbours.local_shifts]  =  Build_Graph_Neighbours_Array_v2(row,col,[1 Neighbours.I_resol(2:end)]); % Define graph of correlations between pixels
    Neighbours.neighb  = Neighbours.neighb(1:end-1,:);
    neighboursSM = Neighbours.neighb(:,1:Neighbours.indGraph(2,2));
    s_filters="";
    for i=1:length(Neighbours.I_resol)-1
        s_filters=strcat(s_filters, int2str(Neighbours.I_resol(i)), ", ");
    end
    s_filters=strcat(s_filters, int2str(Neighbours.I_resol(end)));    
    for k=1:2 % Background
        for i=1:length(PPP)
            for j=1:length(SBR)
                ppp=PPP(i); sbr= SBR(j);
                inFile=sprintf("%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f.mat", ...
                    selectedScene, s_back{k}, K, downSam,ppp, sbr);
                path = strcat(dataDir,'/',inFile);
                if ~isfile(path)
                    warning(strcat("Could not find file: ",inFile));
                    %continue;
                end
                load(strcat(dataDir,'/',inFile),'Y');
                Y=full(Y);
                Y=reshape(Y,row,col,[]);
                fprintf("Processing file %i/%i\n",n,length(PPP)*length(SBR)*2);

                if useTargetDetect
                    [Dep, Refl] = estimateDepthTof(Y, Neighbours.neighb, params, estimateBackground);        
                else
                    [Dep, Refl] = estimateDepthHist(Y, F, Neighbours.neighb, params, estimateBackground);
                end
                Dscales = Dep * params.Tbin *3*10^8/2;

                disp('Part2: Running Weighted median to filter point cloud')
                Dtype = 'Estimates';%'Hist_Back';
                Start=tic;
                [Dest, Rest, ~, D_uncert, R_uncert] = Fct_Denoise_WMF_v2(Dep,Refl,Neighbours, params,F,Dtype, 0,0,params.CvPC, 3, 1);
                toc(Start);        
                dae=sum(sum(abs(Dest-Dref)))/(row*col);
                DAE(i,j,k)=dae;
                IAE(i,j,k)=sum(sum(abs(Rest-Iref)))/(row*col);

                if(1)
                    fprintf("Saving file %i/%i\n",n,length(PPP)*length(SBR)*2);
                    outFile=sprintf("%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f.mat", selectedScene, s_back{k}, K, downSam, ppp, sbr); 
                    I_resol = Neighbours.I_resol;
                    save(strcat(outDir,'/',outFile), 'Dest', 'Dep', 'Rest', 'Refl', 'I_resol');
                end
                if(0)
                    s_title = sprintf("PPP=%i, SBR=%i, %s, filters=[%s], target det=%i, DAE=%f", ppp, sbr, s_back, s_filters, useTargetDetect, DAE);
                    v_axis = [min(Dref(:)), max(Dref(:))];
                    figure; sgtitle(s_title);
                    subplot(2,3,1);imagesc(Dscales);title('Classical Depth');caxis(v_axis);colorbar
                    subplot(2,3,2);imagesc(Dest);title('Proposed Depth');caxis(v_axis);colorbar
                    subplot(2,3,4);imagesc(Refl /max(max(Refl) ) );title('Classical Intensity');colorbar
                    subplot(2,3,5);imagesc(Rest/max(max(Rest)));title('Proposed Intensity');colorbar
                    subplot(2,3,3);imagesc(reshape(D_uncert(:),row,col));title('D uncertainty');caxis([0, 0.02]);colorbar
                    subplot(2,3,6);imagesc(reshape(R_uncert(:),row,col));title('R uncertainty');caxis([0, 10]);colorbar
                end 
                n=n+1;
            end 
        end
    end
end
