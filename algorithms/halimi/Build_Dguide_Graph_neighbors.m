function [ D_up_bar, D_med2, Lvect]   = Build_Dguide_Graph_neighbors(DnMat,r,c,NeighboursSR,Neighbours,ThreshDep,GuideD,ParamPC,Tbin)


LvectSR      = size(NeighboursSR,2);
Lvect        = size(Neighbours,2);
[N,R0]     = size(DnMat);
Dbin = 3*10^8 /2*Tbin;
switch GuideD
    case 0 % Multiscale median
        
        p = 1;D_up_bar = [];I_up_bar = [];
        for ell=1:R0
            for k=1:Lvect
                D_up_bar(p,:) =   DnMat(Neighbours(:,k),ell)  ; % nd SCale  X  N
                p = p+1;
            end
        end
        D_med   = DnMat';
        D_up_bar(D_up_bar<max(Dbin,0.01)) =  NaN;
        indNoDetec =  sum(~isnan(D_up_bar))<=(1*Lvect);
        
        for ell =1:R0
            % Replace bad points by median
            vect        = zeros(Lvect,N);
            DiffNeighb = squeeze( abs((D_up_bar([(ell-1)*Lvect+1:Lvect*ell],:) - D_up_bar([(ell-1)*Lvect+1],:))))<ThreshDep;
            SimPix = sum( DiffNeighb ,1) <= (sqrt(Lvect)+1) ;%|  (D_up_bar([(ell-1)*Lvect+1],:)==0); % L x N
            vect(DiffNeighb) = NaN;
            D_med(ell,SimPix)  = squeeze(nanmedian(vect(:,SimPix).*D_up_bar([(ell-1)*Lvect+1:Lvect*ell],SimPix))); % L x N
            D_med(ell, indNoDetec) = 0;
        end
        D_med(isnan(D_med))=0;
        
        % %         pp=1;
        % %         for ell =1:R0
        % %             % keyboard
        % %             for k=1:LvectSR
        % %                 dd =  D_med(ell,:)';
        % %                 D_med2(pp,:) =   dd(NeighboursSR(:,k))  ;
        % %                 pp = pp+1;
        % %             end
        % %         end
        
        
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
        % %         keyboard
        %         pp=1;
        %         for ell =1:R0
        %             for k=1:LvectSR
        %                 dd =  D_med(ell,:)';
        %                 D_med2(pp,:) =   dd(NeighboursSR(:,k))  ;
        %                 pp = pp+1;
        %             end
        %         end
end

pp=1;
for ell =1:R0
    for k=1:LvectSR
        dd =  D_med(ell,:)';
        D_med2(pp,:) =   dd(NeighboursSR(:,k))  ;
        pp = pp+1;
    end
end

Lvect=LvectSR;
p = 1;D_up_bar = [];
for ell=1:R0
    for k=1:Lvect
        D_up_bar(p,:) =   DnMat(NeighboursSR(:,k),ell)  ; % nd SCale  N
        p = p+1;
    end
end


