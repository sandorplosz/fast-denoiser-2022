clearvars

PPP = [0.1, 1, 5, 10, 50, 100];
SBR  = [0.1, 1, 5, 10, 50, 100];
useTargetDetect=0;
% F_gauss_sig_6=1,  F_real_proc=2
selirf=2;
init

n=1;
for k=1:2
    for i=1:length(PPP)
        for j=1:length(SBR)
            fname=sprintf("Samples_%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f", selectedScene, s_back(k), K, downSam, PPP(i), SBR(j));
            load(strcat(dataDir, '/', fname,".mat"), "Y");
            Y=full(Y');
            %fprintf("%s\n",outFile);
            fprintf("Saving file %i/%i\n", n, 2*length(PPP)*length(SBR))
            fid = fopen(strcat(dataDir,"/bin/",fname,".bin") , 'w'); 
            if fid==-1
                error("Can't open file: %s", fname);
                break
            end
            fwrite(fid , Y , 'uint16');
            fclose(fid);    
            n=n+1;
        end
    end
end