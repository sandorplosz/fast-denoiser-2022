clearvars
useTargetDetect = 1; 
% F_gauss_sig_6=1,  F_real_proc=2
selirf=1;
run ../init.m
PPP = [0.1, 1, 5, 10, 50 100];
SBR  = [0.1, 1, 5, 10, 50 100];

outDir = strcat('../../results/lindell/',irfs(selirf));

%ppp=0.1; sbr=0.1;
DAE = zeros(length(PPP),length(SBR),2);
%IAE = zeros(length(PPP),length(SBR),2);
SumValid = zeros(length(PPP),length(SBR),2);

n=1;
for k=1:2 % Background
    for i=1:length(PPP)
        for j=1:length(SBR)
            ppp=PPP(i); sbr= SBR(j);
            inFile=sprintf("Samples_%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f.mat", ...
                selectedScene, s_back{k}, K, downSam,ppp, sbr);
            path = strcat(outDir,'/',inFile);
            if ~isfile(path)
                warning(strcat("Could not find file: ",inFile));
                %continue;
            end
            load(strcat(outDir,'/',inFile));

            [s1, s2] = size(Dref);
            r1 = 64 - mod(s1,64);
            r2 = 64 - mod(s2,64);
            r1_l = floor(r1/2);
            r2_l = floor(r2/2);
            im2 = im*params.Tbin*3*10^8/2;
            im2 = im2(r1_l:r1_l+s1-1,r2_l:r2_l+s2-1);        
                 
            mask=im2~=0;
            sumValid = sum(mask(:));
            SumValid(i,j,k) = sumValid;
            %dae=sum(sum(abs(Dest(mask)-Dref(mask))))/(sumValid);
            DAE(i,j,k)=sum(abs(im2(:)-Dref(:)))/sumValid;       
            %ISE(i,j,k)=sum(sum(abs(Iest(mask)-IrefGray(mask))));

            if(0)
                fprintf("Saving file %i/%i\n",n,length(PPP)*length(SBR)*2);
                outFile=sprintf("%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f.mat", selectedScene, s_back{k}, K, downSam, ppp, sbr); 
                save(strcat(outDir,'/',outFile), 'Dep', 'Refl', "DAE", "ISE");
            end            
            n=n+1;
        end 
    end
end

save(strcat(outDir,"/res_calculated.mat"), 'DAE', 'SumValid');