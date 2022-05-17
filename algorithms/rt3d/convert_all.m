clearvars;
addpath("..\..\data_generation\")
selirf=2;
useTargetDetect=0;
run ..\..\algorithms\proposed\init.m

PPP = [0.1, 1, 5, 10, 50, 100];
SBR  = [0.1, 1, 5, 10, 50, 100];

if length(h)>K
    h = h(1:K);
else
    d = zeros(K,1);
    d(1:length(h)) = h;
    h = d;
end

n=1;
for k=1:2 % Background
    for i=1:length(PPP)
        for j=1:length(SBR)
            ppp=PPP(i); sbr= SBR(j);
            filename=sprintf('%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f', ...
                    selectedScene, s_back{k}, K, downSam,ppp, sbr);
            path = strcat(dataDir,'/','Samples_', filename, '.mat');
            if ~isfile(path)
                warning(strcat("Could not find file: ",filename));
                %continue;
            end
            load(path,'Y');
            Y=full(Y);
            Y=reshape(Y,row,col,[]);
            outFile=strcat('./real-time-single-photon-lidar/datasets/',filename); 
            fprintf("Converting file %i/%i\n",n,length(PPP)*length(SBR)*2);
            mat_to_RT3D(Y, h, 1, outFile);
            n=n+1;
        end
    end
end

