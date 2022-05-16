function [Dscales, Rscales, Dest, Rest, D_uncert, R_uncert] = detect(Y,F, Neighbours, useTargetDetect, estimateBackground)
    [row, col, T] = size(Y);
    Alpha  = 1;  % non-informative prior
    Beta   = 1;  % non-informative prior 
    NbrePhoton  = 10;       % Approximation with 10 photons (choose <12); 
    limitC      = (nchoosek(NbrePhoton,round(NbrePhoton/2)));
    ProbPrior   = 0.5*ones(row*col,1);
    sigIRF      = 3;
    Attack   = 3;%7  ; % left side of IRF
    trailing = 26;%50;  % right side of IRF
    IRFw        = max([Attack,trailing]); % width of IRF
    Tbin        = 20*10^(-12); % time sample or bin in seconds
    SLight     =  3*10^8 ;
    ThreshDep   = (-9/198*Tbin*10^(12) + 10+18/198)  *SLight/2*Tbin; %1bin for 200ps, 10bins fo
    CvPC        = 0;%0; % 1:conv,   0: operations on PC
    N=row*col;
    

    if useTargetDetect
        disp('    Creation TOF from histogram!)')
        n = 1;
        Nmax = max(max(sum(Y,3)));
        for j=   1:col
            for i=   1:row
                y = squeeze(Y(i,j,:));
                Tof=[];
                ind = find(y>0);
                for k=1:length(ind)
                    Tof = [Tof kron(ind(k),ones(1,y(ind(k))) )];
                end
                Tof    = Tof(randperm(length(Tof)));            
                tof{n} = Tof;
                n      = n+1;
            end
        end
        Tofmat = zeros(N,Nmax);
        nt     = cellfun(@length,tof);
        for n = 1:N
            %ind = randperm(Nmax,nt(n));
            Tofmat(n,1:length(tof{n})) = sort(tof{n});
        end
    end

    Back=1/T*ones(1,T);
    if(estimateBackground)
        disp('    Estimation of Background shape')
        if ~useTargetDetect
            %%% method 1: In case you   have histograms Y and impulse response F
            [~, DD ] = max(F(1:end-1,1:end-1)' * reshape(Y(:,:,2:end),N,T-1)');

        else                   %%% method 2: In case you do not have histograms Y and only have TOF
            DD=zeros(1,N);
            for n= 1:N
                TOFn  = Tofmat(n,:);
                TOFn  = TOFn(TOFn>1);
                [~, ~, DD(n), ~] = ...
                    detect_AH_NonUnifBack_v7_v2(TOFn   ,sigIRF^2,Alpha,Beta,T,NbrePhoton,ProbPrior(n),limitC, 3, 1/T*ones(1,T), 0);%5)
            end
            clear TOFn Det_TD  DiffLogP R_TD 
        end

        %%%%%%%%%%%%%%%%%%%%
        DD      = DD(:);
        Dist    = abs(DD - DD(Neighbours.neighb));
        Dist2   = sort(Dist,2,'ascend');
        mask    = sum(Dist2(:,2:4),2)>50;   %Dist2(:,2)>5;
        Back    = hist(DD.*mask(:),1:T);
        Back(1) = Back(2);
        Back    = movmean(Back,Attack  + trailing);
        Back    = Back/sum(Back);
    end

    if useTargetDetect
        debug=0;
        Start=tic;
        Det_TD=zeros(1,N); DiffLogP=zeros(1,N);
        Dep=zeros(1,N); Refl=zeros(1,N);
        for n= 1:N                    
            TOFn  = Tofmat(n,:);%(Iter-2)*Group+1:(Iter-1)*Group);
            TOFn  = TOFn(TOFn>1);
            TOFn=sort(TOFn,'descend');
            [Det_TD(n), DiffLogP(n), Dep(n),Refl(n)] = ...
                detect_AH_NonUnifBack_v7_v2(TOFn   ,sigIRF^2,Alpha,Beta,T,NbrePhoton,ProbPrior(n),limitC, 3, Back,debug);%5)
        end
        toc(Start);
        if(1)
            dep=reshape(Dep,row,col);
            figure; imagesc(dep);
        end
        indd=(sum(Y(:,:,2:end),3) == 1);  % AH-mod
        Det_TD(indd) = 1;            % AH-mod    
        Dep  = Dep(:).*Det_TD(:);
        ind     = (Det_TD(:)==0);
        r0      = reshape(Y(:,:,1),N,1);
        %R_TD(ind)    = r0(ind);
        Refl(ind)    = 10-r0(ind); Refl(r0>=10)=1; % AH-mod    
    else
        % Running simple Cross Correlation
        % v1
        Yt = reshape(Y,[],T);
        [~,d1]=max((Yt*F),[],2);

        % v2
        Yt = reshape(Y,[],T)';
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
        Yf = fft(Yt);
        XcorrMatrix = ifft(Yf.*hf,'symmetric'); 
        if 1
            XcorrMatrix = XcorrMatrix - Back';        
        elseif(estimateBackground)        
            XcorrMatrix=XcorrMatrix';
            Ysort     = sort(XcorrMatrix); % N  T
            ProfileT  = median(Ysort(1:floor(0.2*N),:,:)); % 1 T L            
            ProfileT  = (ProfileT-mean(ProfileT,2)); % 1 T L
            ProfileN  = median(XcorrMatrix(:,1:T),2); %N 1 L
            XcorrMatrix  = max(XcorrMatrix -  max(0,ProfileT + ProfileN), 0)';
        end      
        [~, Dep] = max(XcorrMatrix);     
        %[a(:,1), a(:,2)] = sort(XcorrMatrix(:,90),'descend');
        Refl=zeros(row,col);
        n=1;
        for k=1:col
            for l=1:row
                Refl(l,k) = sum(Yt(max(1,Dep(n)-7):min(Dep(n)+50,T),n)); %7 is left width of IRF, 50 is the right width
                n = n +1;
            end
        end   
        if(0)
            ind=Dep(:)<10;
            Dep(ind)=0;
            r0 = reshape(Y(:,:,1),row*col,1);
            Refl(ind) = 10-r0(ind); Refl(r0>=10)=1; % AH-mod
        end
    end

    Dep=reshape(Dep,row,col);
    Dscales = Dep *Tbin*3*10^8/2;% AH-mod
    Rscales = reshape(Refl,row,col);
    %figure; imagesc(Dscales);

    disp('Part2: Running Weighted median to filter point cloud')
    Dtype = 'Estimates';%'Hist_Back';
    Start=tic;
    [Dest, Rest, ~, D_uncert, R_uncert] = Fct_Denoise_WMF_v2(Dep,Refl,Neighbours.I_resol,F,IRFw,Tbin,Attack  ,trailing,Neighbours.neighb,Neighbours.indGraph,Neighbours.local_shifts,Dtype,ThreshDep,0,0,CvPC);
    toc(Start);

end