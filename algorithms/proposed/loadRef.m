inFile=sprintf("%s_%s_K_%i_DownS_%i*.mat", ...
            selectedScene, backName, K, downSam);
files=dir(strcat(dataDir, '/', inFile));
load(strcat(files(1).folder,'/',files(1).name));
