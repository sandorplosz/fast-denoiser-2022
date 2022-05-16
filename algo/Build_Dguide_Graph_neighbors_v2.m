function [ D_med, Lvect, D_up_bar]   = Build_Dguide_Graph_neighbors(DnMat,r,c,Neighbours,local_shifts,ThreshDep,GuideD,ParamPC,Tbin,CvPC)


Lvect        = size(Neighbours,2);
[N,R0]     = size(DnMat);
Dbin = 3*10^8 /2*Tbin; 

        if(CvPC==1) % Do convolutions
            diamK = sqrt(Lvect);radiusK = floor(diamK/2);
            pos = zeros(diamK);
            for i=1:Lvect %v=pos;v(i)=1; 
                posM(radiusK+1+local_shifts(i,1),radiusK+1+local_shifts(i,2),i) = 1;  
            end
            Dmat = zeros(r+2*radiusK,c+2*radiusK,R0);
            Dmat(radiusK+1:end-radiusK,radiusK+1:end-radiusK,:) = reshape(DnMat,[r,c,R0]);
            
            for ell=1:R0
                D_up_bar22(:,:,(ell-1)*Lvect+1:ell*Lvect) = convn(posM,Dmat(:,:,ell));
            end
            D_up_bar22 = D_up_bar22(2*radiusK+1:end-2*radiusK, 2*radiusK+1:end-2*radiusK, :);
            D_up_bar = reshape(D_up_bar22,N,[])';
            
        else % Do point wise operations
            p = 1;D_up_bar = [];I_up_bar = [];pp=1;
            for ell=1:R0
                for k=1:Lvect
                    D_up_bar(p,:) =   DnMat(Neighbours(:,k),ell)  ; % nd SCale  X  N
                    p = p+1;
                end
            end
        end
        
        switch GuideD
    case 0 % Multiscale median
        D_med   = DnMat';
        D_up_bar(D_up_bar<max(Dbin,0.01)) =  NaN;
        indNoDetec =  sum(~isnan(D_up_bar))<=(1*Lvect);
        
        for ell =1:R0
            % Replace bad points by median
            vect        = ones(Lvect,N);
            % Shouldn't this be the opposite? e.g. >ThreshDep
            DiffNeighb = squeeze( abs((D_up_bar([(ell-1)*Lvect+1:Lvect*ell],:) - D_up_bar([(ell-1)*Lvect+1],:))))<ThreshDep;
            SimPix = sum( DiffNeighb ,1) <= (sqrt(Lvect)+1) ;%|  (D_up_bar([(ell-1)*Lvect+1],:)==0); % L x N
            vect(DiffNeighb) = NaN;
            D_med(ell,SimPix)  = squeeze(nanmedian(vect(:,SimPix).*D_up_bar([(ell-1)*Lvect+1:Lvect*ell],SimPix))); % L x N
            D_med(ell, indNoDetec) = 0;
        end
        D_med(isnan(D_med))=0;
        
        
    case 1 % PC denoise
        DD = DnMat;
        DD(DnMat<max(Dbin,0.01)) =  NaN;
        PC         = [repmat(kron((1:c)', ones(r,1)),3,1), (DD(:))/Dbin  , repmat(kron(ones(c,1),(r:-1:1)'),3,1)];
        PC         = pointCloud(PC);
        if(ParamPC==0)
            [PCest,In,Out] = pcdenoise(PC,'threshold',0.01,'NumNeighbors',15);
        else
            [PCest,In,Out] = pcdenoise(PC,'threshold',0.05,'NumNeighbors',5);
        end
        D_med      = DnMat;
        D_med(DnMat<max(Dbin,0.01)) = 0;
        D_med(Out) = 0;
        D_med      = D_med';
end
D_up_bar(isnan(D_up_bar)) = 0;


