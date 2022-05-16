% Load lindell
clearvars
useTargetDetect = 0;
% F_gauss_sig_6=1,  F_real_proc=2
selirf=2;
init
PPP = [0.1, 1, 5, 10, 30];
SBR  = [0.1, 1, 5, 10, 30];

path_lindell = strcat('/home/ps2014/Development/matlab/fast_denoiser_2022/results_lindell/',irfs(selirf),'/');
path_abde_tci = strcat('/home/ps2014/Development/matlab/fast_denoiser_2022/results_abde_tci/',irfs(selirf),'/');
path_prop = strcat('/home/ps2014/Development/matlab/fast_denoiser_2022/results_algo/',irfs(selirf),'/');
DAE_lindell = zeros(length(PPP),length(SBR),2);
IAE_lindell = zeros(length(PPP),length(SBR),2);
DAE_tci = zeros(length(PPP),length(SBR),2);
DAE_prop = zeros(length(PPP),length(SBR),2);
DAE_xcorr = zeros(length(PPP),length(SBR),2);
%ppp=0.1; sbr=5;

for k=1:2
    for i=1:length(PPP)
        for j=1:length(SBR)
            ppp=PPP(i); sbr= SBR(j);
            inFile=sprintf("%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f", ...
                selectedScene, s_back{k}, K, downSam,ppp, sbr);

            fpath = strcat(path_lindell,'Samples_',inFile,'.mat');
            if ~isfile(fpath)
                warning(strcat("Could not find file1: ",inFile));            
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

            fpath = strcat(path_abde_tci,'/',inFile,'.mat');
            if ~isfile(fpath)
                warning(strcat("Could not find file2: ",inFile));            
            else
                res_tci = load(fpath);
                DAE_tci(i,j,k)=sum(sum(abs(res_tci.Dest-Dref)))/(row*col);
            end

            fpath = strcat(path_prop,'/',inFile,'.mat');
            if ~isfile(fpath)
                warning(strcat("Could not find file3: ",inFile));            
            else
                res_prop = load(fpath);
                res_prop.Dep2=res_prop.Dep*params.Tbin *3*10^8/2;
                DAE_prop(i,j,k)=sum(sum(abs(res_prop.Dest-Dref)))/(row*col);
                DAE_xcorr(i,j,k)=sum(sum(abs(res_prop.Dep2-Dref)))/(row*col);
            end
        end
    end
end

s_filters="";
for i=1:length(res_prop.I_resol)-1
    s_filters=strcat(s_filters, int2str(res_prop.I_resol(i)), ", ");
end
s_filters=strcat(s_filters, int2str(res_prop.I_resol(end)));

s_back = ["Uniform Background", "Gamma Background"];

%res1=load('/home/ps2014/Development/matlab/fast_denoiser_2022/results_algo/error_wmf_dguide.mat');
%res2=load('/home/ps2014/Development/matlab/fast_denoiser_2022/results_algo/error_maxfilter_noguide_mean2.mat');
temp = [DAE_xcorr(:); DAE_lindell(:); DAE_prop(:); DAE_tci(:)];
temp=log10(temp);
c_axis = [min(temp), max(temp)];
c_axis=[-2.8, 0.4];

figure;
n=1;
for k=1:2
    subplot(2,4,n);
    contourf(log10(PPP),log10(SBR),log10(DAE_xcorr(:,:,k)'));
    title("Classic XCorr")
    ylabel({s_back{k},'SBR'});xlabel('Photons');caxis(c_axis)
    xticks([-1 0 1 2]);xticklabels({'10^{-1}','1','10^1','10^2'})
    yticks([-1 0 1 2]);yticklabels({'10^{-1}','1','10^1','10^2'})
    colorbar
    n=n+1;    
    subplot(2,4,n);
    contourf(log10(PPP),log10(SBR),log10(DAE_tci(:,:,k)'));
    title("Halimi")
    ylabel('SBR');xlabel('Photons');caxis(c_axis)
    xticks([-1 0 1 2]);xticklabels({'10^{-1}','1','10^1','10^2'})
    yticks([-1 0 1 2]);yticklabels({'10^{-1}','1','10^1','10^2'})
    colorbar
    n=n+1;
    subplot(2,4,n);
    contourf(log10(PPP),log10(SBR),log10(DAE_lindell(:,:,k)'));
    title("Lindell")
    ylabel('SBR');xlabel('Photons');caxis(c_axis)
    xticks([-1 0 1 2]);xticklabels({'10^{-1}','1','10^1','10^2'})
    yticks([-1 0 1 2]);yticklabels({'10^{-1}','1','10^1','10^2'})
    colorbar
    n=n+1;
    subplot(2,4,n);
    contourf(log10(PPP),log10(SBR),log10(DAE_prop(:,:,k)'));
    title('Proposed');
    ylabel('SBR');xlabel('Photons');caxis(c_axis)
    xticks([-1 0 1 2]);xticklabels({'10^{-1}','1','10^1','10^2'})
    yticks([-1 0 1 2]);yticklabels({'10^{-1}','1','10^1','10^2'})
    colorbar
    n=n+1;
end

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
res_tci = load(strcat('/home/ps2014/Development/matlab/fast_denoiser_2022/results_abde_tci/',inFile));

%v_axis = [min(Dref(:)) max(Dref(:))];
v_axis = [0.9*min(Dref(:)), 1.1*max(Dref(:))];
%figure; imagesc(Dref); caxis(v_axis);
%figure; imagesc(res_lindell.im2); caxis(v_axis); colorbar
props = {'ML',0.06, 'MR', 0.03, 'MT', 0.1, 'MB', 0.07};

% Replace subaxis with subfigure if not available
s_title = sprintf("PPP=%i, SBR=%i, %s, filters=[%s], target det=%i", ppp, sbr, s_back, s_filters, useTargetDetect);
figure; sgtitle(s_title);
subaxis(1,4,1,props{:}); imagesc(Dref); caxis(v_axis); colorbar; title('Reference Depth');
subaxis(1,4,2); imagesc(res_tci.Dest); caxis(v_axis); colorbar; title('TCI')
subaxis(1,4,3); imagesc(res_lindell.im2); caxis(v_axis); colorbar; title('Lindell')
subaxis(1,4,4); imagesc(Dest); caxis(v_axis); colorbar; title('Proposed')
end
