function [Dep, Refl] = estimateDepthTof(Y, neighbours, p, estimateBackground, dopostproc)

    if exist('dopostproc','var') ~=1
        dopostproc = true;
    end    

    [row, col, T] = size(Y);
    N=row*col;
   
    disp('    Creation TOF from histogram')
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
    for n = 1:N
        %ind = randperm(Nmax,nt(n));
        Tofmat(n,1:length(tof{n})) = sort(tof{n});
    end
    
    if(estimateBackground)
        disp('    Estimation of Background shape')
               %%% method 2: In case you do not have histograms Y and only have TOF
        DD=zeros(1,N);
        Start=tic;
        for n= 1:N
            TOFn  = Tofmat(n,:);
            TOFn  = TOFn(TOFn>1);
            [~, ~, DD(n), ~] = ...
                detect_AH_NonUnifBack_v7_v2(TOFn   ,p.sigIRF^2,p.Alpha,p.Beta,T,p.NbrePhoton,p.ProbPrior(n),p.limitC, 3, 1/T*ones(1,T), 0);%5)
        end
        toc(Start);
        clear TOFn     
        %dd=reshape(DD,row,col); figure;imagesc(dd); figure;imagesc(td_dep1)
        %run mres_td.m
        %sum(DD~=td_dep1(:)')
        %DD=td_dep1(:)';

        %%%%%%%%%%%%%%%%%%%%
        DD      = DD(:);
        Dist    = abs(DD - DD(neighbours));
        Dist2   = sort(Dist,2,'ascend');
        mask    = sum(Dist2(:,2:4),2)>50;   %Dist2(:,2)>5;
        Back    = hist(DD.*mask(:),1:T);
        Back(1) = Back(2);
        a=Back;
        %Back    = movmean(Back,[1,2]);
        Back    = movmean(Back,[p.Attack, p.trailing]);               
        Back    = Back/sum(Back);
        %temp=[a;Back;td_backg']; 
        %sum(abs(Back-td_backg')>0.001)
        %Back=td_backg';
        
    else
        Back=1/T*ones(1,T);
    end
    
    debug=0;    
    Det_TD=zeros(1,N); DiffLogP=zeros(1,N);
    Dep=zeros(1,N); Refl=zeros(1,N);
    disp('    Running detection')
    Start=tic;
    for n= 1:N                    
        TOFn  = Tofmat(n,:);%(Iter-2)*Group+1:(Iter-1)*Group);
        TOFn  = TOFn(TOFn>1);
        TOFn=sort(TOFn,'descend');
        [Det_TD(n), DiffLogP(n), Dep(n),Refl(n)] = ...
            detect_AH_NonUnifBack_v7_v2(TOFn   ,p.sigIRF^2,p.Alpha,p.Beta,T,p.NbrePhoton,p.ProbPrior(n),p.limitC, 3, Back,debug);%5)
    end
    toc(Start);
    if(0)
        dep=reshape(Dep,row,col);
        figure; imagesc(dep);
        [dv, di]=sort(abs(dep(:)-td_depth(:)),'descend');
        n=di(5);
        [dep(n),td_depth(n)] 
        sum(dep(:)~=td_dep2(:))
        figure;imagesc(abs(dep-td_dep2))
    end

    if(dopostproc)
        indd=(sum(Y(:,:,2:end),3) == 1);  % AH-mod
        Det_TD(indd) = 1;            % AH-mod    
        Dep  = Dep(:).*Det_TD(:);
        ind     = (Det_TD(:)==0);
        r0      = reshape(Y(:,:,1),N,1);
        %R_TD(ind)    = r0(ind);
        Refl(ind)    = 10-r0(ind); Refl(r0>=10)=1; % AH-mod    
    end
    
    Dep=reshape(Dep,row,col);
    %Dep = Dep *Tbin*3*10^8/2;% AH-mod
    Refl = reshape(Refl,row,col);    
    
end