% Based on Abdi's script Build_Synth_Data_TCI2021_robustMS_UnifGammaB.m

function [ Dref, IrefGray] = build_synth(D_HR, I_HR, F, K, PPP, SBR, choix, downSam, outFile)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Median filter: Reconstruct missing parts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Row, Col]    = size(D_HR);
    [indI, indJ] = find(D_HR==0);
    DiffD = (min(cat(3,abs(D_HR - circshift(D_HR,[-1,0])),abs(D_HR - circshift(D_HR,[1,0]))...
        ,abs(D_HR - circshift(D_HR,[0,1])) , abs(D_HR - circshift(D_HR,[0,-1]))),[],3));
    [indI1, indJ1] = find(DiffD > 50);
    [indI] = [indI ; indI1];
    [indJ] = [indJ ; indJ1];

    step         = 10;
    D0_HR        = D_HR;
    for n=1:length(indI)
        xi  = max(1,indI(n)-step):min(Row,indI(n)+step);
        yj  = max(1,indJ(n)-step):min(Col,indJ(n)+step);
        vvv = D_HR(xi,yj);
        D0_HR(indI(n),indJ(n)) = round(median(vvv(:)));
    end
    %figure; subplot(2,1,1);imagesc(D_HR); subplot(2,1,2);imagesc(D0_HR)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Downsampling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dref     =  max(1,D0_HR(2:downSam:end,1:downSam:end));
    Iref     =  max(1,I_HR(2:downSam:end,1:downSam:end,:));

    [Row, Col]  = size(Dref);
    N = Row*Col;

    Iref = Iref/max(Iref(:));
    IrefGray = rgb2gray(Iref);  
    IrefGray = IrefGray/mean(IrefGray(:));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 3D data cube: 3 wavelengths
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Sg = zeros(Row,Col,K);
    Dref = K-Dref-100; %% !!! Shift histogram to the left !!!
    for i=1:Row            
        for j=1:Col
            Sg(i,j,Dref(i,j)) = (IrefGray(i,j));
        end
    end

    %figure; I1=IrefGray(:,:,1);
    %scatter3(kron([1:Col]',ones(Row,1)),Dref(:),kron(ones(Col,1),[Row:-1:1]'),1,I1(:),'.'),
    
    Sg_Conv =  reshape((F*reshape(Sg(:,:,:),[N K])')', [Row Col K]) ;   
    Sg_Poiss = zeros(Row,Col,K);
    Yg_Poiss = zeros(Row,Col,K);    
    n=1;

    switch choix
        case 'GENERATE_NOISY_CUBES'
            for LevP = 1:length(PPP)                
                for LevSBR = 1:length(SBR)
                    Lev_S = SBR(LevSBR) .*PPP(LevP) ./(1+SBR(LevSBR));
                    Lev_B = PPP(LevP)- Lev_S;     
                    fprintf("Sampling case PPP=%f,SBR=%f\n",PPP(LevP),SBR(LevSBR));
                    Sg_Poiss = poissrnd(Sg_Conv*Lev_S);                    
                    Yg_Poiss = Sg_Poiss + poissrnd(Lev_B*ones(Row, Col, K)/K);  
                    fileName=sprintf(outFile, PPP(LevP), SBR(LevSBR));
                    Y=sparse(reshape(Yg_Poiss,N,[]));
                    fprintf("Saving file %i/%i\n",n,length(PPP)*length(SBR));
                    save(fileName, 'Y');     
                    n=n+1;
                end
            end
            
        case 'GENERATE_NOISY_CUBES_GAMMA'
            %%Shape B
            t = 1:1500; shapeB=circshift(0.2 * t.^1.5 .* exp(-0.01*t),[0 -20]); 
            shapeB = shapeB(1:K);
            shapeB = shapeB / sum(shapeB) * K;
            %figure;plot(shapeB)
            %saveas(gcf,'gamma_back_shape.fig');
            %shapeB  = gampdf(1:K,2,30)* K; 
            shapeB  = repmat(shapeB(:)',N,1);
            shapeB  = reshape(shapeB,[Row,Col,K]);
            for LevP = 1:length(PPP)                
                for LevSBR = 1:length(SBR)
                    Lev_S = SBR(LevSBR) .*PPP(LevP) ./(1+SBR(LevSBR));
                    Lev_B = PPP(LevP)- Lev_S;
                    fprintf("Sampling case PPP=%f,SBR=%f\n",PPP(LevP),SBR(LevSBR));
                    Sg_Poiss = poissrnd(Sg_Conv*Lev_S);
                    Yg_Poiss = Sg_Poiss + poissrnd(Lev_B*shapeB/K);
                    fileName=sprintf(outFile, PPP(LevP), SBR(LevSBR));
                    Y=sparse(reshape(Yg_Poiss,N,[]));
                    fprintf("Saving file %i/%i\n",n,length(PPP)*length(SBR));
                    save(fileName, 'Y');     
                    n=n+1;
                end
            end    
        otherwise
            fprintf("Invalid option: %s", choix );
    end
end
