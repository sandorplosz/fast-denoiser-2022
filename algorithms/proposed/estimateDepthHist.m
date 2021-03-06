function [Dep, Refl] = estimateDepthHist(Y, F, neighbours, p, estimateBackground, useCenterOfMass, Dref) 

    if exist('useCenterOfMass','var') ~=1
        useCenterOfMass = false;
    end

    [row, col, T] = size(Y);
    N=row*col;    
    
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

    %%%%TEMP CODE
    %h=single(h);
    %Yt=single(Yt);
    %%%

    %h=log(h+eps); % the log is optional
    hf = fft(h);
    Yf = fft(Yt);
    XcorrMatrix = ifft(Yf.*hf,'symmetric'); 
    %[xv, xi] = sort(XcorrMatrix(:,1156),'descend');
    %         
    if(estimateBackground==1)        
        [~, DD] = max(XcorrMatrix);    
        DD      = DD(:);
        %fprintf("DAE before: %f\n", log10(sum(abs(DD(:)*p.Tbin*3*10^8/2-Dref(:)))/N));
        Dist    = abs(DD - DD(neighbours));
        Dist2   = sort(Dist,2,'ascend');
        mask    = sum(Dist2(:,2:8),2)>50;   %Dist2(:,2)>5;
        Back    = hist(DD.*mask(:),1:T);
        Back(1) = Back(2);
        Back    = movmean(Back,p.Attack  + p.trailing);
        Back    = Back/sum(Back);  
        %XcorrMatrix = XcorrMatrix - Back'; 
        XcorrMatrix = max(XcorrMatrix - Back',0); 
        %[~,a]=max(XcorrMatrix- Back'); a=a*p.Tbin*3*10^8/2;
        %[~,b]=max(max(XcorrMatrix- Back',2)); b=b*p.Tbin*3*10^8/2;
        %log10(sum(abs(a(:)*p.Tbin*3*10^8/2-Dref(:)))/N)
        %log10(sum(abs(b(:)*p.Tbin*3*10^8/2-Dref(:)))/N)
    elseif(estimateBackground==2)        
        XcorrMatrix=XcorrMatrix';
        Ysort     = sort(XcorrMatrix); % N  T
        ProfileT  = median(Ysort(1:floor(0.2*N),:,:)); %          
        ProfileT  = (ProfileT-mean(ProfileT,2)); 
        ProfileN  = median(XcorrMatrix(:,1:T),2);
        XcorrMatrix  = max(XcorrMatrix -  max(0,ProfileT + ProfileN), 0)';       
    elseif(estimateBackground==3)    
        Back=zeros(T,1);
        count=zeros(T,1);
        [M, DD] = max(XcorrMatrix);    
        DD      = DD(:);
        Dist    = abs(DD - DD(neighbours));
        Dist2   = sort(Dist,2,'ascend');
        mask    = sum(Dist2(:,2:4),2)>50;   %Dist2(:,2)>5;
        for i=1:N
            if mask(i)
                a=DD(i);
                Back(a)=Back(a)+M(i);
                count(a)=count(a)+1;
            end
        end
        Back=Back./count;
        Back(isnan(Back))=0;
        Back    = movmean(Back,p.Attack  + p.trailing);
        %Back    = Back/sum(Back);  
        XcorrMatrix = XcorrMatrix - Back;  
    end  

    Refl=zeros(row,col);    

    if(~useCenterOfMass)
        [~, Dep] = max(XcorrMatrix);
        n=1;
        for k=1:col
            for l=1:row
                Refl(l,k) = sum(Yt(max(1,Dep(n)-p.Attack):min(Dep(n)+p.trailing,T),n));
                n = n +1;
            end
        end           
    else
        Dep=zeros(N, 1);
        [~, xmax] = max(XcorrMatrix);    
        for i=1:N
            v=max(xmax(i)-p.Attack,1):min(xmax(i)+p.trailing,T);
            Dep(i)=v*Yt(v,i)/sum(eps+Yt(v,i));  
            Refl(i)=sum(Yt(v,i));
        end
    end
    %Dscales = Dep*p.Tbin *3*10^8/2;
    %log10(sum(abs(Dref(:)-Dscales(:)))/N)
    %[a(:,1), a(:,2)] = sort(XcorrMatrix(:,90),'descend');
  
    %ind=Dep(:)<=1;
    %Dep(ind)=0;  
    if(0)
        r0 = reshape(Y(:,:,1),row*col,1);
        Refl(ind) = 10-r0(ind); Refl(r0>=10)=1; % AH-mod
    end
   
    Dep=reshape(Dep,row,col);
    Refl = reshape(Refl,row,col);    
end
