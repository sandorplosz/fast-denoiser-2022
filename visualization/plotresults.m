clearvars
useTargetDetect = 0;
% F_gauss_sig_6=1,  F_real_proc=2
selirf=2;
init
PPP = [0.1, 1, 5, 10, 50, 100];
SBR  = [0.1, 1, 5, 10, 50, 100];

path_class = strcat('../results/classic/',irfs(selirf),'/');
path_lindell = strcat('../results/lindell/',irfs(selirf),'/');
path_halimi = strcat('../results/halimi/',irfs(selirf),'/');
path_rt3d = strcat('../results/rt3d/',irfs(selirf),'/');
path_prop_mat = strcat('../results/proposed_matlab/',irfs(selirf),'/');
path_prop = strcat('../results/proposed/',irfs(selirf),'/');

DAE_class = zeros(length(PPP),length(SBR),2);
DAE_lindell = zeros(length(PPP),length(SBR),2);
IAE_lindell = zeros(length(PPP),length(SBR),2);
DAE_halimi = zeros(length(PPP),length(SBR),2);
DAE_rt3d = zeros(length(PPP),length(SBR),2);
DAE_prop = zeros(length(PPP),length(SBR),2);
DAE_prop_mat = zeros(length(PPP),length(SBR),2);

rt3d_valid = zeros(length(PPP),length(SBR),2);

warning('on')

%load(strcat('../results/classic/', irfs(selirf),'/res_calculated.mat'));

for k=1:2
    for i=1:length(PPP)
        for j=1:length(SBR)
            ppp=PPP(i); sbr= SBR(j);
            inFile=sprintf("%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f", ...
                selectedScene, s_back{k}, K, downSam,ppp, sbr);

            fpath = strcat(path_class,'/',inFile,'.mat');
            if ~isfile(fpath)
                warning(strcat("Could not find file: ",fpath));            
            else
                res_class = load(fpath);
                DAE_class(i,j,k)=sum(sum(abs(res_class.Dep-Dref_orig)))/(row*col)*params.Tbin *3*10^8/2;
            end                      

            if(1)
                fpath = strcat(path_lindell,'Samples_',inFile,'.mat');
                if ~isfile(fpath)
                    warning(strcat("Could not find file: ",fpath));            
                else
                    res_lindell = load(fpath);
                    [s1, s2] = size(Dref);
                    r1 = 64 - mod(s1,64);
                    r2 = 64 - mod(s2,64);
                    r1_l = floor(r1/2);
                    r2_l = floor(r2/2);
                    res_lindell.im2 = res_lindell.im*params.Tbin*3*10^8/2;
                    res_lindell.im2=res_lindell.im2(r1_l:r1_l+s1-1,r2_l:r2_l+s2-1);        
                    DAE_lindell(i,j,k)=sum(sum(abs(res_lindell.im2-Dref)))/(row*col);
                    %IAE_lindell(i,j)=sum(sum(abs(Rest-Iref)))/(row*col);
                end
            end

            fpath = strcat(path_halimi,'/',inFile,'.mat');
            if ~isfile(fpath)
                warning(strcat("Could not find file: ",fpath));            
            else
                res_halimi = load(fpath);
                DAE_halimi(i,j,k)=sum(sum(abs(res_halimi.Dest-Dref)))/(row*col);
            end

            fpath = strcat(path_rt3d,'/',inFile,'.mat');
            if ~isfile(fpath)
                warning(strcat("Could not find file: ",fpath));            
            else
                res_rt3d = load(fpath);                
                DAE_rt3d(i,j,k)=res_rt3d.dae*params.Tbin*3*10^8/2;
                rt3d_valid(i,j,k) = res_rt3d.sumValid/(row*col);
            end            

            fpath = strcat(path_prop_mat,'/',inFile,'.mat');
            if ~isfile(fpath)
                warning(strcat("Could not find file: ",fpath));            
            else
                res_prop_mat = load(fpath);
                DAE_prop_mat(i,j,k)=sum(sum(abs(res_prop_mat.Dest-Dref)))/(row*col);
            end            
    
            inFile2=sprintf("%s_app.m", strrep(inFile,'.','_'));
            fpath = strcat(path_prop,inFile2);
            if ~isfile(fpath)
                warning(strcat("Could not find file: ", fpath));            
            else
                run(fpath)
                DAE_prop(i,j,k)=sum(sum(abs(den_depth-Dref)))/(row*col);
                %DAE_xcorr(i,j,k)=sum(sum(abs(res_prop.Dep2-Dref)))/(row*col);
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

%res1=load('/home/ps2014/Development/matlab/fast_denoiser_2022/results_algo/error_wmf_dguide.mat');
%res2=load('/home/ps2014/Development/matlab/fast_denoiser_2022/results_algo/error_maxfilter_noguide_mean2.mat');
temp = [DAE_class(:); DAE_lindell(:); DAE_prop(:); DAE_halimi(:)];
temp=log10(temp);
temp=temp(~isinf(temp));
c_axis = [min(temp), max(temp)];
%c_axis=[-2.8, 0.4];

figure;
n=1;
for k=1:2
    subplot(2,6,n);
    contourf(log10(PPP),log10(SBR),log10(DAE_class(:,:,k)'));
    title("Classic XCorr")
    ylabel({s_back_full{k},'SBR'});xlabel('Photons');caxis(c_axis)
    xticks([-1 0 1 2]);xticklabels({'10^{-1}','1','10^1','10^2'})
    yticks([-1 0 1 2]);yticklabels({'10^{-1}','1','10^1','10^2'})
    colorbar
    n=n+1;    

    subplot(2,6,n);
    contourf(log10(PPP),log10(SBR),log10(DAE_halimi(:,:,k)'));
    title("Halimi")
    ylabel('SBR');xlabel('Photons');caxis(c_axis)
    xticks([-1 0 1 2]);xticklabels({'10^{-1}','1','10^1','10^2'})
    yticks([-1 0 1 2]);yticklabels({'10^{-1}','1','10^1','10^2'})
    colorbar
    n=n+1;

    subplot(2,6,n);
    contourf(log10(PPP),log10(SBR),log10(DAE_lindell(:,:,k)'));
    title("Lindell")
    ylabel('SBR');xlabel('Photons');caxis(c_axis)
    xticks([-1 0 1 2]);xticklabels({'10^{-1}','1','10^1','10^2'})
    yticks([-1 0 1 2]);yticklabels({'10^{-1}','1','10^1','10^2'})
    colorbar
    n=n+1;

    subplot(2,6,n);
    contourf(log10(PPP),log10(SBR),log10(DAE_rt3d(:,:,k)'));
    title("RT3D")
    ylabel('SBR');xlabel('Photons');caxis(c_axis)
    xticks([-1 0 1 2]);xticklabels({'10^{-1}','1','10^1','10^2'})
    yticks([-1 0 1 2]);yticklabels({'10^{-1}','1','10^1','10^2'})
    colorbar
    n=n+1;    

    subplot(2,6,n);
    contourf(log10(PPP),log10(SBR),log10(DAE_prop_mat(:,:,k)'));
    title('Proposed (matlab)');
    ylabel('SBR');xlabel('Photons');caxis(c_axis)
    xticks([-1 0 1 2]);xticklabels({'10^{-1}','1','10^1','10^2'})
    yticks([-1 0 1 2]);yticklabels({'10^{-1}','1','10^1','10^2'})
    colorbar   
    n=n+1;

    subplot(2,6,n);
    contourf(log10(PPP),log10(SBR),log10(DAE_prop(:,:,k)'));
    title('Proposed');
    ylabel('SBR');xlabel('Photons');caxis(c_axis)
    xticks([-1 0 1 2]);xticklabels({'10^{-1}','1','10^1','10^2'})
    yticks([-1 0 1 2]);yticklabels({'10^{-1}','1','10^1','10^2'})
    colorbar
    n=n+1;
end

figure
contourf(log10(PPP),log10(SBR),rt3d_valid(:,:,k)');
title('RT3D proportion of valid pixels');
ylabel('SBR');xlabel('Photons');caxis([0,1])
xticks([-1 0 1 2]);xticklabels({'10^{-1}','1','10^1','10^2'})
yticks([-1 0 1 2]);yticklabels({'10^{-1}','1','10^1','10^2'})
colorbar

if(0)
    figure;
    n=1;
    for k=1:2
        subplot(2,1,n);
        contourf(log10(PPP),log10(SBR),log10(DAE_lindell(:,:,k)));
        title("Lindell's")
        ylabel('SBR');xlabel('Photons');caxis([-3 0.5])
        xticks([-1 0 1]);xticklabels({'10^{-1}','10^0','10^1'})
        yticks([-1 0 1 2]);yticklabels({'10^{-1}','1','10^1','10^2'})
        colorbar
        n=n+1;
    end
end

if(0)
%Load TCI
res_halimi = load(strcat('/home/ps2014/Development/matlab/fast_denoiser_2022/results_abde_tci/',inFile));

%v_axis = [min(Dref(:)) max(Dref(:))];
v_axis = [0.9*min(Dref(:)), 1.1*max(Dref(:))];
%figure; imagesc(Dref); caxis(v_axis);
%figure; imagesc(res_lindell.im2); caxis(v_axis); colorbar
props = {'ML',0.06, 'MR', 0.03, 'MT', 0.1, 'MB', 0.07};

% Replace subaxis with subfigure if not available
s_title = sprintf("PPP=%i, SBR=%i, %s, filters=[%s], target det=%i", ppp, sbr, s_back, s_filters, useTargetDetect);
figure; sgtitle(s_title);
subaxis(1,4,1,props{:}); imagesc(Dref); caxis(v_axis); colorbar; title('Reference Depth');
subaxis(1,4,2); imagesc(res_halimi.Dest); caxis(v_axis); colorbar; title('TCI')
subaxis(1,4,3); imagesc(res_lindell.im2); caxis(v_axis); colorbar; title('Lindell')
subaxis(1,4,4); imagesc(Dest); caxis(v_axis); colorbar; title('Proposed')
end
