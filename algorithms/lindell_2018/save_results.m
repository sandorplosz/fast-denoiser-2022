clearvars
useTargetDetect = 1; 
% F_gauss_sig_6=1,  F_real_proc=2
selirf=1;
init
PPP = [0.1, 1, 5, 10, 50 100];
SBR  = [0.1, 1, 5, 10, 50 100];

outDir = strcat('../../results/lindell/',irfs(selirf));

%ppp=0.1; sbr=0.1;
DSE = zeros(length(PPP),length(SBR),2);
ISE = zeros(length(PPP),length(SBR),2);
SumValid = zeros(length(PPP),length(SBR),2);

n=1;
for k=1:2 % Background
    for i=1:length(PPP)
        for j=1:length(SBR)
            ppp=PPP(i); sbr= SBR(j);
            inFile=sprintf("%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f.mat", ...
                selectedScene, s_back{k}, K, downSam,ppp, sbr);
            path = strcat(outDir,'/',inFile);
            if ~isfile(path)
                warning(strcat("Could not find file: ",inFile));
                %continue;
            end
            load(strcat(outDir,'/',inFile));
         
            mask=Dest~=0;
            SumValid(i,j,k) = sum(mask(:));
            %dae=sum(sum(abs(Dest(mask)-Dref(mask))))/(sumValid);
            DSE(i,j,k)=sum(sum(abs(Dest(mask)-Dref(mask))));
            ISE(i,j,k)=sum(sum(abs(Iest(mask)-IrefGray(mask))));

            if(0)
                fprintf("Saving file %i/%i\n",n,length(PPP)*length(SBR)*2);
                outFile=sprintf("%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f.mat", selectedScene, s_back{k}, K, downSam, ppp, sbr); 
                save(strcat(outDir,'/',outFile), 'Dep', 'Refl', "DAE", "ISE");
            end            
            n=n+1;
        end 
    end
end

save(strcat(outDir,"/res_calculated.mat"), 'DSE', 'ISE', 'SumValid');