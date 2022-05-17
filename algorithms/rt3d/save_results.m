ppp=10; sbr=10; 
k=1;

filename=sprintf('%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f', ...
        selectedScene, s_back{k}, K, downSam,ppp, sbr);
path = ['./real-time-single-photon-lidar/output_' filename '/frame0_w0.ply'];

if ~isfile(path)
    error("Could not find file: %s",filename);
    %continue;
end

ptCloud = pcread(path);
%pcshow(ptCloud);
points = ptCloud.Location;
Dep=zeros(row,col);
Refl=zeros(row,col);
%ind=(points(:,2)-1)*row+points(:,3);
%dep(ind)=points(:,1);
for i = 1:size(points,1)
    r=row-points(i,3)+1;
    c=col-points(i,2)+1;
    p = points(i,1);
    if(Dep(r,c)==0 || abs(Dref_orig(r,c)-p)<abs(Dref_orig(r,c)-Dep(r,c)))    
        Dep(r,c) = p;
        Refl(r,c)=ptCloud.Intensity(i);
    end
end
mean(abs(Dep(:)-Dref_orig(:)))

ca=[min(Dref_orig(:)), max(Dref_orig(:))];
figure;imagesc(Dep); caxis(ca);
figure;imagesc(Dref_orig);
figure;imagesc(Refl);
figure;imagesc(IrefGray);


% Scale_rat=0.01 numPPP=1 neighb=1 -- NO POINTS
% Scale_rat=0.01 numPPP=2 neighb=1 -- NO POINTS
% Scale_rat=0.01 numPPP=1 neighb=2 -- NO POINTS
% Scale_rat=0.01 numPPP=2 neighb=2 -- NO POINTS
% Scale_rat=1 numPPP=1 neighb=1 -- 22.88
% Scale_rat=1 numPPP=2 neighb=1 -- 767.2741
% Scale_rat=1 numPPP=1 neighb=2 -- 20.6076

save(['..\matlab\fast_denoiser_2022\results_rt3d\' filename '.mat'], "Dep")
