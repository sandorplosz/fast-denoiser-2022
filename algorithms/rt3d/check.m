PPP = [0.1, 1, 5, 10, 50, 100];
SBR  = [0.1, 1, 5, 10, 50, 100];

n = 0;
for k=1:2 % Background
    for i=1:length(PPP)
        for j=1:length(SBR)
            ppp=PPP(i); sbr= SBR(j);
            filename=sprintf('%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f', ...
                    selectedScene, s_back{k}, K, downSam,ppp, sbr);
            file = strcat('./real-time-single-photon-lidar/output_', filename, '/frame0_w0.ply');
            if ~isfile(file)
                %warning(strcat("Could not find file: ",filename));
                n=n+1;
                %continue;
            end
        end
    end
end

fprintf("%i/%i files not found\n", n, length(PPP)*length(SBR)*2);