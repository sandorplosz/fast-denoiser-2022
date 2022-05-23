clearvars
clc
scenedir="./raw";
scenes = { 'Reindeer', 'Art', 'Plastic', 'Moebius', 'Laundry', ...
            'Dolls', 'Bowling1', 'Books' };
selectedScene=2;
downSam = 2;
T = 1024;

D_HR   = double(imread(strcat(scenedir, '/', scenes{selectedScene}, '/disp1.png')));
I_HR   = double(imread(strcat(scenedir, '/', scenes{selectedScene}, '/view1.png')));

PPP = [0.1, 1, 5, 10, 50 100];
SBR  = [0.1, 1, 5, 10, 50 100];
Lev_S = SBR .*PPP ./(1+SBR);
Lev_B = PPP - Lev_S;

backChoices = ["GENERATE_NOISY_CUBES", "GENERATE_NOISY_CUBES_GAMMA"];
backNames = ["UnifBack", "GammaBack"];

irf=2; % 1: Real, 2: Gaussian

if(irf==1)
    irfname="F_real_proc";
else
    sigma2=6^2;
    irfname=sprintf("F_gauss_sig_%i",sqrt(sigma2));
end

if(1) % Generate IRF
    if(irf==1)
        %% True IRF
        load F_real2_100s
        F = processF(F,T);
        h=F(:,ceil(size(F,2)/2));
        F_window=getIRFWindow(h,0.1);
        save(strcat(irfname,'.mat'),'F', 'F_window');
        fid = fopen(sprintf('./%s_%i_%i.bin',irfname,F_window(1),F_window(2)) , 'w');
        fwrite(fid , h , 'float32');
        fclose(fid);      
    else    
        h= 1/sqrt(2*pi*sigma2)* exp(- ((1:T)-round(T/2)).^2/2/(sigma2));
        h = h(:)/sum(h);
        [~,attack] = max(h);
        l = length(h);
        F=zeros(l,l);    
        for i=1:l
            F(:,i)  = circshift(h,[i-attack 0]);
        end
        F_window=getIRFWindow(h,0.1);
        F = F(1:T,1:T);
        save(strcat(irfname,'.mat'),'F', 'F_window');
        fid = fopen(strcat('./', irfname, '.bin') , 'w');
        fwrite(fid , h , 'float32');
        fclose(fid);    
    end
else
    load(strcat(irfname,'.mat'));
end

for i=1:2
    outFile=sprintf("%s/Samples_%s_%s_K_%i_DownS_%i_PPP_%%.3f_SBR_%%.3f.mat", irfname, scenes{selectedScene}, backNames(i), T, downSam);
    [Dref, IrefGray] = build_synth(D_HR, I_HR, F, T, PPP, SBR, backChoices(i), downSam, outFile);
end

save(strcat(scenes{selectedScene}, '_ref_img.mat'), 'Dref', 'IrefGray');

if(0)    
    h = F(:,round(T/2));   % T'x1 Select one IRF or maybe you have it in a file
    %%% Create a vector h of size Tx1, with max at 0
    if length(h)>size(Yt,1)
        h = h(1:size(Yt,1));
    else
        d = zeros(size(Yt,1),1);
        d(1:length(h)) = h;
        h = d;
    end
    %h=h_orig;
    h = h/sum(h); % Tx1
    [~,attack] = max(h);
    h = circshift(h,-attack); 
    %%%
    h  = flipud(h); % flip h in preparation to FFT fft(log(h+eps))    
    %h=log(h+eps); % the log is optional
    hf = fft(h);
    sg=Sg(:,:,:,2,2);
    Yt = reshape(sg,[],T)';
    Yf = fft(Yt);
    XcorrMatrix = ifft(Yf.*hf,'symmetric'); 
    [~, DD] = max(XcorrMatrix);    
    DD=reshape(DD,row,col);
    sum(DD(:)~=Dref(:))
    dae=sum(abs(DD(:)-Dref_orig(:)))/N;
    a=[DD(:),Dref(:)];
end

%[R,C] = size(Dref);
%N=R*C;

if(0)
    n=1;
    for i=1:length(PPP)
        for j=1:length(SBR)
            fprintf("Saving file %i/%i\n",n,length(PPP)*length(SBR));
            n=n+1;
            outFile=sprintf("Samples_%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f.mat", scenes{selectedScene}, backName, T, downSam, PPP(i), SBR(j));
            Y=Ys(:,:,:,i,j);
            Y=sparse(reshape(Y,N,[]));
            %fprintf("%s\n",outFile);
            save(outFile, 'Dref', 'Y', 'Iref');
        end
    end
end







