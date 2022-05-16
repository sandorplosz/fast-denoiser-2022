clearvars

PPP = [0.1, 1, 5, 10, 30];
SBR  = [0.1, 1, 5, 10, 30];
useTargetDetect=0;
% F_gauss_sig_6=1,  F_real_proc=2
selirf=2;
init

n=1;
for k=1:2
    for i=1:length(PPP)
        for j=1:length(SBR)
            fname=sprintf("Samples_%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f", selectedScene, s_back(k), K, downSam, PPP(i), SBR(j));
            load(strcat(irfs(selirf),'/',fname,".mat"), "Y");
            Y=full(Y');
            %fprintf("%s\n",outFile);
            fprintf("Saving file %i/%i\n", n, 2*length(PPP)*length(SBR))
            fid = fopen(strcat(irfs(selirf),"/bin/",fname,".bin") , 'w');            
            fwrite(fid , Y , 'uint16');
            fclose(fid);    
            n=n+1;
        end
    end
end