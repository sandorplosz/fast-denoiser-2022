clearvars;
addpath("..\..\data_generation\")
selirf=2;
useTargetDetect=0;
run ..\..\algorithms\proposed\init.m

ppp=10; sbr=10; 
k=1;

if length(h)>K
    h = h(1:K);
else
    d = zeros(K,1);
    d(1:length(h)) = h;
    h = d;
end

filename=sprintf('%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f', ...
        selectedScene, s_back{k}, K, downSam,ppp, sbr);
path = strcat(dataDir,'/','Samples_', filename, '.mat');
if ~isfile(path)
    warning(strcat("Could not find file: ",filename));
    %continue;
end
load("../../data_generation/Art_ref_img.mat");
load(path,'Y');
Y=full(Y);
Y=reshape(Y,row,col,[]);

outFile=strcat('./real-time-single-photon-lidar/datasets/',filename); 

%% convert it
%mat_to_RT3D(filename);
mat_to_RT3D(Y, h, 0.01, outFile);
