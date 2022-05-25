clearvars
% F_gauss_sig_6=1,  F_real_proc=2
selirf=1;
run ../init.m
PPP = [0.1, 1, 5, 10, 50 100];
SBR  = [0.1, 1, 5, 10, 50 100];

outDir = strcat('../../results/classic/',irfs(selirf));

%ppp=0.1; sbr=0.1;
DAE = zeros(length(PPP),length(SBR),2);
IAE = zeros(length(PPP),length(SBR),2);

h = h/sum(h);
[~,attack] = max(h);
h = circshift(h,-attack); 
h  = flipud(h); % flip h in preparation to FFT fft(log(h+eps))  
hf = fft(h);

n=1;
for k=1:2 % Background
    for i=1:length(PPP)
        for j=1:length(SBR)
            ppp=PPP(i); sbr= SBR(j);
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
            Y=full(Y');            
            fprintf("Processing file %i/%i\n",n,length(PPP)*length(SBR)*2);
            Yf = fft(Y);
            
            XcorrMatrix = ifft(Yf.*hf,'symmetric'); 
            [~, Dep] = max(XcorrMatrix);
            Refl=sum(Y);
            if(0)
                dep=reshape(Dep,row,col);
                ref=reshape(Refl,row,col);
                figure;imagesc(dep);
                figure;imagesc(ref);
            end

            Dep=Dep*params.Tbin*3*10^8/2;
         
            dae=sum(abs(Dep(:)-Dref(:)))/(row*col);
            DAE(i,j,k)=dae;
            IAE(i,j,k)=sum(abs(Refl(:)-Iref(:)))/mean(Iref(:));

            if(0)
                fprintf("Saving file %i/%i\n",n,length(PPP)*length(SBR)*2);
                outFile=sprintf("%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f.mat", selectedScene, s_back{k}, K, downSam, ppp, sbr); 
                save(strcat(outDir,'/',outFile), 'Dep', 'Refl', "DAE", "IAE");
            end            
            n=n+1;
        end 
    end
end

save(strcat(outDir,"/res_calculated.mat"), 'DAE', 'IAE');