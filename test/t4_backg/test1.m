clearvars
useTargetDetect = 0;
init;
T=K;
Dtype = 'Estimates';%'Hist_Back';    

outDir = '/home/ps2014/Development/matlab/fast_denoiser_2022/results_algo';
addpath('/home/ps2014/Development/matlab/fast_denoiser_2022/algo')
load /home/ps2014/Development/matlab/fast_denoiser_2022/data_generation/F_real2_100s
F = processF(F,K);
h = F(:,round(T/2));   % T'x1 Select one IRF or maybe you have it in a file
h = h/sum(h); % Tx1
[~,attack] = max(h);
h = circshift(h,-attack); 
%%%
h  = flipud(h); % flip h in preparation to FFT fft(log(h+eps))    
%h=log(h+eps); % the log is optional
hf = fft(h);

Neighbours.I_resol = [1 5 7 9] ;  % size of spatial correlations Requires
%Neighbours.I_resol = [1 7 13] ;  % size of spatial correlations Requires
[Neighbours.neighb, Neighbours.indGraph, Neighbours.local_shifts]  =  Build_Graph_Neighbours_Array_v2(row,col,[1 Neighbours.I_resol(2:end)]); % Define graph of correlations between pixels
Neighbours.neighb  = Neighbours.neighb(1:end-1,:);
neighboursSM = Neighbours.neighb(:,1:Neighbours.indGraph(2,2));
s_filters="";
for i=1:length(Neighbours.I_resol)-1
    s_filters=strcat(s_filters, int2str(Neighbours.I_resol(i)), ", ");
end
s_filters=strcat(s_filters, int2str(Neighbours.I_resol(end)));

ppp=1;
sbr=100;
back = 1;
%estimateBack = 3;

inFile=sprintf("%s_%s_K_%i_DownS_%i_PPP_%.3f_SBR_%.3f.mat", ...
    selectedScene, s_back{back}, K, downSam,ppp, sbr);
path = strcat(dataDir,'/',inFile);
if ~isfile(path)
    warning(strcat("Could not find file: ",inFile));
    %continue;
end
load(strcat(dataDir,'/',inFile),'Y');
Y=full(Y);
Y=reshape(Y,row,col,[]);
Ysum=squeeze(sum(sum(Y,1),2));
figure;plot(Ysum);

[row, col, T] = size(Y);
N=row*col;    

Yt = reshape(Y,[],T)';
Yf = fft(Yt);
XcorrMatrix = ifft(Yf.*hf,'symmetric'); 
[~, DD] = max(XcorrMatrix);    
DD      = DD(:);
fprintf("DAE before: %f\n", log10(sum(abs(DD*params.Tbin*3*10^8/2-Dref(:)))/N));
for i=2:length(Neighbours.I_resol)
    Dist    = abs(DD - DD(Neighbours.neighb(:,1:Neighbours.indGraph(2,i))));
    Dist2   = sort(Dist,2,'ascend');
    mask    = sum(Dist2(:,2:4),2)>50;   %Dist2(:,2)>5;
    Back    = hist(DD.*mask(:),1:T);
    Back(1) = Back(2);
    Back    = movmean(Back,params.Attack  + params.trailing);
    Back    = Back/sum(Back);  
    figure;plot(Back);title(sprintf("%i",Neighbours.indGraph(2,i)))
    %XCM = max(0, XcorrMatrix-Back');
    XCM = XcorrMatrix-Back';
    [~, Dep] = max(XCM);   
    Refl=zeros(row,col);
    n=1;
    for k=1:col
        for l=1:row
            Refl(l,k) = sum(Yt(max(1,Dep(n)-7):min(Dep(n)+50,T),n));
            n = n +1;
        end
    end   
    ind=Dep(:)<=1;   
    Dep(ind)=0;  
    Dep=reshape(Dep,row,col);
    Refl = reshape(Refl,row,col);    
    Dscales = Dep * params.Tbin *3*10^8/2;
    fprintf("DAE with filter %i after backgr.: %f", Neighbours.I_resol(i), log10(sum(abs(Dscales(:)-Dref(:)))/N));
    [Dest, Rest, ~, ~, ~] = Fct_Denoise_WMF_v3(Dep,Refl,Neighbours.I_resol,F,params.IRFw, params.Tbin, ...
            params.Attack, params.trailing, Neighbours.neighb, Neighbours.indGraph, Neighbours.local_shifts,Dtype, params.ThreshDep,0,0,params.CvPC, 0, 1);
    dae=log10(sum(sum(abs(Dest-Dref)))/(row*col));    
    fprintf(", after denoise: %i\n",dae);
end

[Dep, Refl] = estimateDepthHist(Y, F, Neighbours.neighb, params, 2,Dref);
Dscales = Dep * params.Tbin *3*10^8/2;
dae=log10(sum(sum(abs(Dscales-Dref)))/(row*col));
fprintf("DAE after est: %f\n", dae);

     
