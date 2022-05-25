clearvars
useTargetDetect = 0;
% F_gauss_sig_6=1,  F_real_proc=2
selirf=1;
run ../algorithms/init.m
PPP = [0.1, 1, 5, 10, 50, 100];
SBR  = [0.1, 1, 5, 10, 50, 100];
show_prop_matlab=0;

path_class = strcat('../results/classic/',irfs(selirf),'/');
path_lindell = strcat('../results/lindell/',irfs(selirf),'/res_calculated.mat');
path_halimi = strcat('../results/halimi/',irfs(selirf),'/');
path_rt3d = strcat('../results/rt3d/',irfs(selirf),'/results_processed.mat');
path_prop_mat = strcat('../results/proposed_matlab/',irfs(selirf),'/');
path_prop = strcat('../results/proposed/',irfs(selirf),'/');

DAE_class = zeros(length(PPP),length(SBR),2);
IAE_class = zeros(length(PPP),length(SBR),2);
DAE_halimi = zeros(length(PPP),length(SBR),2);
IAE_halimi = zeros(length(PPP),length(SBR),2);
DAE_prop = zeros(length(PPP),length(SBR),2);
IAE_prop = zeros(length(PPP),length(SBR),2);
DAE_prop_mat = zeros(length(PPP),length(SBR),2);

%res_classic = load(path_class);
%res_classic.DAE=res_classic.DAE*params.Tbin*3*10^8/2;
res_lindell = load(path_lindell);
res_rt3d = load(path_rt3d);
res_rt3d.DAE=res_rt3d.DAE*params.Tbin*3*10^8/2;
res_rt3d.SumValid=res_rt3d.SumValid/(row*col);
res_rt3d.DAE(res_rt3d.SumValid<0.2)=0;

warning('on')

%load(strcat('../results/classic/', irfs(selirf),'/res_calculated.mat'));

for k=1:2
    for i=1:length(PPP)
        for j=1:length(SBR)
            ppp=PPP(i); sbr= SBR(j);
            Lev_S = sbr*ppp/(1+sbr);
            Iref = IrefGray * Lev_S;

            inFile=sprintf("%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f", ...
                selectedScene, s_back{k}, K, downSam,ppp, sbr);

            fpath = strcat(path_class,'/',inFile,'.mat');
            if ~isfile(fpath)
                warning(strcat("Could not find file: ",fpath));            
            else
                res_class = load(fpath);
                res_class.Dep=res_class.Dep*params.Tbin *3*10^8/2;
                DAE_class(i,j,k)=sum(abs(res_class.Dep(:)-Dref(:)))/(row*col);
                IAE_class(i,j,k)=sum(abs(res_class.Refl(:)-Iref(:)))/mean(Iref(:));
            end                      
        
%             fpath = strcat(path_lindell,'Samples_',inFile,'.mat');
%             if ~isfile(fpath)
%                 warning(strcat("Could not find file: ",fpath));            
%             else
%                 res_lindell = load(fpath);
%                 [s1, s2] = size(Dref);
%                 r1 = 64 - mod(s1,64);
%                 r2 = 64 - mod(s2,64);
%                 r1_l = floor(r1/2);
%                 r2_l = floor(r2/2);
%                 res_lindell.im2 = res_lindell.im*params.Tbin*3*10^8/2;
%                 res_lindell.im2=res_lindell.im2(r1_l:r1_l+s1-1,r2_l:r2_l+s2-1);        
%                 DAE_lindell(i,j,k)=sum(sum(abs(res_lindell.im2-Dref)))/(row*col);
%                 %IAE_lindell(i,j)=sum(sum(abs(Rest-Iref)))/(row*col);
%             end
            
            fpath = strcat(path_halimi,'/',inFile,'.mat');
            if ~isfile(fpath)
                warning(strcat("Could not find file: ",fpath));            
            else
                res_halimi = load(fpath);
                DAE_halimi(i,j,k)=sum(abs(res_halimi.Dest(:)-Dref(:)))/(row*col);
                IAE_halimi(i,j,k)=sum(abs(res_halimi.Iest(:)-Iref(:)))/mean(Iref(:));
            end

%             if(0)
%                 fpath = strcat(path_rt3d,'/',inFile,'.mat');
%                 if ~isfile(fpath)
%                     warning(strcat("Could not find file: ",fpath));            
%                 else
%                     res_rt3d = load(fpath);                
%                     DAE_rt3d(i,j,k)=res_rt3d.dae*params.Tbin*3*10^8/2;
%                     rt3d_valid(i,j,k) = res_rt3d.sumValid/(row*col);
%                 end            
%             end

            if(show_prop_matlab)
                fpath = strcat(path_prop_mat,'/',inFile,'.mat');
                if ~isfile(fpath)
                    warning(strcat("Could not find file: ",fpath));            
                else
                    res_prop_mat = load(fpath);
                    DAE_prop_mat(i,j,k)=sum(sum(abs(res_prop_mat.Dest-Dref)))/(row*col);
                end            
            end
    
            inFile2=sprintf("%s_app.m", strrep(inFile,'.','_'));
            fpath = strcat(path_prop,inFile2);
            if ~isfile(fpath)
                warning(strcat("Could not find file: ", fpath));            
            else
                run(fpath)
                DAE_prop(i,j,k)=sum(abs(den_depth(:)-Dref(:)))/(row*col);
                IAE_prop(i,j,k)=sum(abs(den_refl(:)-Iref(:)))/mean(Iref(:));
            end  
        end
    end
end

% s_filters="";
% for i=1:length(res_prop.I_resol)-1
%     s_filters=strcat(s_filters, int2str(res_prop.I_resol(i)), ", ");
% end
% s_filters=strcat(s_filters, int2str(res_prop.I_resol(end)));
s_filters="1 3 7";
s_back_full = ["Uniform Background", "Gamma Background"];

dae_scale = [DAE_class(:); res_lindell.DAE(:); DAE_prop(:); DAE_halimi(:); res_rt3d.DAE(:)];
dae_scale=log10(dae_scale);
dae_scale=dae_scale(~isinf(dae_scale));
c_axis = [min(dae_scale), max(dae_scale)];
%c_axis=[-2.8, 0.4];

numfig = 5 + 1*(show_prop_matlab~=0);
ticksx = [-1 0 1 2]; labelx = {'10^{-1}','1','10^1','10^2'};
ticksy = [-1 0 1 2]; labely = {'10^{-1}','1','10^1','10^2'};

figure;
tl=tiledlayout(2,numfig, "TileSpacing","compact", "Padding","compact");
xlabel(tl,'Average number of photons per pixel')
ylabel(tl,'SBR')
for k=1:2
    nexttile
    contourf(log10(PPP),log10(SBR),log10(DAE_class(:,:,k)'));
    if(k==1), title("Classic XCorr"); end
    ylabel(s_back_full{k},'FontWeight','Bold');
    caxis(c_axis)
    xticks(ticksx);xticklabels(labelx)
    yticks(ticksy);yticklabels(labely)   

    nexttile
    contourf(log10(PPP),log10(SBR),log10(DAE_halimi(:,:,k)'));
    if(k==1), title("Halimi"); end
    caxis(c_axis)
    xticks(ticksx);xticklabels(labelx)
    yticks(ticksy);yticklabels(labely)

    nexttile
    contourf(log10(PPP),log10(SBR),log10(res_lindell.DAE(:,:,k)'));
    if(k==1), title("Lindell"); end
    caxis(c_axis)
    xticks(ticksx);xticklabels(labelx)
    yticks(ticksy);yticklabels(labely)

    nexttile
    contourf(log10(PPP),log10(SBR),log10(res_rt3d.DAE(:,:,k)'));
    if(k==1), title("RT3D"); end
    caxis(c_axis)
    xticks(ticksx);xticklabels(labelx)
    yticks(ticksy);yticklabels(labely)

    if show_prop_matlab
        nexttile
        contourf(log10(PPP),log10(SBR),log10(DAE_prop_mat(:,:,k)'));
        if(k==1), title('Proposed (matlab)'); end
        caxis(c_axis)
        xticks(ticksx);xticklabels(labelx)
        yticks(ticksy);yticklabels(labely)
        %colorbar       
    end

    nexttile
    contourf(log10(PPP),log10(SBR),log10(DAE_prop(:,:,k)'));
    if(k==1), title('Proposed'); end
    caxis(c_axis)
    xticks(ticksx);xticklabels(labelx)
    yticks(ticksy);yticklabels(labely)
    
end
cb = colorbar;
cb.Layout.Tile = 'east';

if(0)
    figure
    tl=tiledlayout(2,1, "TileSpacing","compact");
    ylabel(tl,'SBR');xlabel(tl,'Photons per pixel');
    title(tl,'RT3D proportion of returned pixels');
    for k=1:2
        nexttile
        contourf(log10(PPP),log10(SBR),res_rt3d.SumValid(:,:,k)');        
        caxis([0,1])
    xticks(ticksx);xticklabels(labelx)
    yticks(ticksy);yticklabels(labely)
    end
    cb = colorbar;
    cb.Layout.Tile = 'east';
end


iae_scale = [IAE_class(:); IAE_prop(:); IAE_halimi(:)];
iae_scale=log10(iae_scale);
iae_scale=iae_scale(~isinf(iae_scale));

c_axis = [min(iae_scale), max(iae_scale)];
figure;
tl=tiledlayout(2,3, "TileSpacing","compact", "Padding","compact");
xlabel(tl,'Average number of photons per pixel')
ylabel(tl,'SBR')
for k=1:2
    nexttile
    contourf(log10(PPP),log10(SBR),log10(IAE_class(:,:,k)'));
    if(k==1), title("Classic XCorr"); end
    ylabel(s_back_full{k},'FontWeight','Bold');
    caxis(c_axis)
    xticks(ticksx);xticklabels(labelx)
    yticks(ticksy);yticklabels(labely)   

    nexttile
    contourf(log10(PPP),log10(SBR),log10(IAE_halimi(:,:,k)'));
    if(k==1), title("Halimi"); end
    caxis(c_axis)
    xticks(ticksx);xticklabels(labelx)
    yticks(ticksy);yticklabels(labely)

    nexttile
    contourf(log10(PPP),log10(SBR),log10(IAE_prop(:,:,k)'));
    if(k==1), title("Proposed"); end
    caxis(c_axis)
    xticks(ticksx);xticklabels(labelx)
    yticks(ticksy);yticklabels(labely)    
end
