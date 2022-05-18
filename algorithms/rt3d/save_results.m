selirf=2;
useTargetDetect=0;
run ../proposed/init.m

PPP = [0.1, 1, 5, 10, 50, 100];
SBR  = [0.1, 1, 5, 10, 50, 100];

n=1;
for k=1:2 % Background
    for i=1:length(PPP)
        for j=1:length(SBR)
            ppp=PPP(i); sbr= SBR(j);

            filename=sprintf('%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f', ...
                    selectedScene, s_back{k}, K, downSam,ppp, sbr);
            path = ['../../results/rt3d/pc/output_' filename '/frame0_w0.ply'];
            
            if ~isfile(path)
                error("Could not find file: %s",filename);
                continue;
            end

            fprintf("Processing file %i/%i\n",n,length(PPP)*length(SBR)*2);
            
            ptCloud = pcread(path);
            %pcshow(ptCloud);
            points = ptCloud.Location;
            Dep=zeros(row,col);
            Refl=zeros(row,col);
            for l = 1:size(points,1)
                r=row-points(l,3)+1;
                c=col-points(l,2)+1;
                p = points(l,1);
                %If there are multiple points per pixel, we store the closest one to
                %the reference!
                if(Dep(r,c)==0 || abs(Dref_orig(r,c)-p)<abs(Dref_orig(r,c)-Dep(r,c)))    
                    Dep(r,c) = p;
                    Refl(r,c)=ptCloud.Intensity(l);
                end
            end
            dae=mean(abs(Dep(:)-Dref_orig(:)));
            iae=mean(abs(Refl(:)-IrefGray(:)))/mean(IrefGray(:));
            
            if(0)
                ca=[min(Dref_orig(:)), max(Dref_orig(:))];
                figure;imagesc(Dep); caxis(ca);
                figure;imagesc(Dref_orig);
                figure;imagesc(Refl);
                figure;imagesc(IrefGray);
            end
            
            %save(['../../results/rt3d/' filename '.mat'], 'Dep', 'Refl', 'dae', 'iae');
            save(['../../results/rt3d/' filename '.mat'], 'dae', 'iae');
            n=n+1;
        end
    end
end
