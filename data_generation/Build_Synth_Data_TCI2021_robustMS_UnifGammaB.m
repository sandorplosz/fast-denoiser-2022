%% Generate synthetic data - TCI Rev1
%% this is the data used for the TCI paper after revision
clear all
clc
choice = 0; %0  one peak,   1: multipeaks
scenedir="C:\LocalStore\ps2014\Development\lindell_2018_code\middlebury\raw";
scenes = { 'Reindeer', 'Art', 'Plastic', 'Moebius', 'Laundry', ...
            'Dolls', 'Bowling1', 'Books' };
selectedScene=1;
downSam = 2;
ppplength = 5;
K = 300;

outFile=sprintf("%s_GT_3Wavelengths_1peak_185_232_%i_IndD.mat", scenes{selectedScene}, K);
D_HR   = double(imread(strcat(scenedir, '\', scenes{selectedScene}, '\disp1.png')));
I_HR   = double(imread(strcat(scenedir, '\', scenes{selectedScene}, '\view1.png')));

Choice{1}  =   'GENERATE_CLEAN_CUBES';%   
Choice{2}  =   'GENERATE_NOISY_CUBES';%      Unif Back
Choice{3}  =   'GENERATE_NOISY_CUBES_GAMMA';%     Gamma back 


for Indchoix = 3
    choix = Choice{Indchoix};
    switch choix
        case 'GENERATE_CLEAN_CUBES'
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
            c_axis = [min(D_HR(:)) max(D_HR(:))];
            figure
            subplot(2,1,1);imagesc(D_HR);caxis(c_axis);colorbar
            subplot(2,1,2);imagesc(D0_HR);caxis(c_axis);colorbar
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Downsampling
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Dref     =  max(1,D0_HR(2:downSam:end,1:downSam:end));
            Iref     =  max(1,I_HR(2:downSam:end,1:downSam:end,:));
            
            [Row Col]  = size(Dref);
            %            Iref = Iref/sum(Iref(:));
            Iref = Iref/max(Iref(:));
            
            IrefGray = rgb2gray(Iref);  
            Iref = Iref/mean(Iref(:));
%             for i=1:3  Iref(:,:,i)  = Iref(:,:,i)/mean(IrefGray(:)) ; end
            IrefGray = IrefGray/mean(IrefGray(:));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% 3D data cube: 3 wavelengths
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            S1 = zeros(Row,Col,K);
            S2 = zeros(Row,Col,K);
            S3 = zeros(Row,Col,K);
            Sg = zeros(Row,Col,K);
            Dref = K-Dref;

            for i=1:Row
                i
                for j=1:Col
                    S1(i,j,Dref(i,j)) = (Iref(i,j,1));
                    S2(i,j,Dref(i,j)) = (Iref(i,j,2));
                    S3(i,j,Dref(i,j)) = (Iref(i,j,3));
                    Sg(i,j,Dref(i,j)) = (IrefGray(i,j));
                end
            end
            figure;imagesc(sum(S1,3))
            
            figure; I1=IrefGray(:,:,1);
            scatter3(kron([1:Col]',ones(Row,1)),Dref(:),kron(ones(Col,1),[Row:-1:1]'),1,I1(:),'.'),
            [Row Col K] = size(S1);
            
%             figure;imagesc(sum(S1~=0,3));
            %
            %  figure; I1=I(:,:,1);
            %  scatter3(kron([1:Col]',ones(Row,1)),Dref(:),kron(ones(Col,1),[Row:-1:1]'),1,I1(:),'.'),
            
%             save Art_GT_3Wavelengths_1peak_185_232_300.mat S1 S2 S3 Sg Dref Iref IrefGray            
            save(outFile, 'S1', 'S2', 'S3', 'Sg', 'Dref', 'Iref', 'IrefGray');
            
        case 'GENERATE_NOISY_CUBES'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%  Noise
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% Generate noisy data/different Counts per pixel
%             load Art_GT_3Wavelengths_1peak_185_232_300.mat S1 S2 S3 Sg Dref Iref IrefGray
            load(outFile, 'S1', 'S2', 'S3', 'Sg', 'Dref', 'Iref', 'IrefGray');

            [Row Col K] = size(S1);
            N           = Row*Col;
            
            %% True IRF
            load F_real2_100s
            F = F(:,100); F(F<(0.01*max(F(:))))=0;
            h = [zeros(100,1); F(2:2:end); zeros(100,1)];
            h  = h/sum(h);
            [vvv, MaxF] = max(h); clear F
            for i=1:length(h)
                F(:,i)  = circshift(h,[i-MaxF 0]);
            end
            F = F(1:K,1:K);
            
            %% SBR, PPP levels
            pppY = logspace(log10(0.1),log10(1000),ppplength);
            SBR  = logspace(log10(0.01),log10(100),ppplength);
            
            %% Convolution in Times by IRF
            S1_Conv =  reshape((F*reshape(S1(:,:,:),[N K])')', [Row Col K]) ;
            S2_Conv =  reshape((F*reshape(S2(:,:,:),[N K])')', [Row Col K]) ;
            S3_Conv =  reshape((F*reshape(S3(:,:,:),[N K])')', [Row Col K]) ;
            Sg_Conv =  reshape((F*reshape(Sg(:,:,:),[N K])')', [Row Col K]) ;
            
            
            
            for LevP = 1:length(pppY)
                LevP
                for LevSBR = 1:length(SBR)
                    Lev_S = SBR(LevSBR) .*pppY(LevP) ./(1+SBR(LevSBR));
                    Lev_B = pppY(LevP)- Lev_S;
                    SBRmat(LevP,LevSBR)  = SBR(LevSBR);
                    pppYmat(LevP,LevSBR) = pppY(LevP);
                    
                    %                     S1_GT(:,:,:,LevP,LevSBR)   = S1*Lev_S;
                    %                     S2_GT(:,:,:,LevP,LevSBR)   = S2*Lev_S;
                    %                     S3_GT(:,:,:,LevP,LevSBR)   = S3*Lev_S;
                    %                     Sg_GT(:,:,:,LevP,LevSBR)   = Sg*Lev_S;
                    
                    S1_Poiss(:,:,:,LevP,LevSBR) = poissrnd(S1_Conv*Lev_S);
                    S2_Poiss(:,:,:,LevP,LevSBR) = poissrnd(S2_Conv*Lev_S);
                    S3_Poiss(:,:,:,LevP,LevSBR) = poissrnd(S3_Conv*Lev_S);
                    Sg_Poiss(:,:,:,LevP,LevSBR) = poissrnd(Sg_Conv*Lev_S);
                    
                    Y1_Poiss(:,:,:,LevP,LevSBR) =  S1_Poiss(:,:,:,LevP,LevSBR) + poissrnd(Lev_B*ones(Row, Col, K)/K);
                    Y2_Poiss(:,:,:,LevP,LevSBR) =  S2_Poiss(:,:,:,LevP,LevSBR) + poissrnd(Lev_B*ones(Row, Col, K)/K);
                    Y3_Poiss(:,:,:,LevP,LevSBR) =  S3_Poiss(:,:,:,LevP,LevSBR) + poissrnd(Lev_B*ones(Row, Col, K)/K);
                    Yg_Poiss(:,:,:,LevP,LevSBR) =  Sg_Poiss(:,:,:,LevP,LevSBR) + poissrnd(Lev_B*ones(Row, Col, K)/K);
                    
                end
            end
%             Name = ['Art_Noisy_3Wav_1peak_185_232_300_5PPPy_5SBR_UnifBack.mat'];
              Name = ['Art_Noisy_3Wav_1peak_185_232_300_5PPPy_5SBR_UnifBack_IndD.mat'];

            save(Name,'SBRmat','pppYmat','pppY','SBR', 'F',  'Dref', 'Iref', 'IrefGray',...
                'S1','S2','S3','Sg', ...
                'Y1_Poiss','Y2_Poiss','Y3_Poiss','Yg_Poiss',...
                'S1_Poiss','S2_Poiss','S3_Poiss','Sg_Poiss','-v7.3')
             
            
            
            
        case 'GENERATE_NOISY_CUBES_GAMMA'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%  Noise
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%  Noise
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% Generate noisy data/different Counts per pixel
%             load Art_GT_3Wavelengths_1peak_185_232_300.mat S1 S2 S3 Sg Dref Iref IrefGray
            load Art_GT_3Wavelengths_1peak_185_232_300_IndD.mat S1 S2 S3 Sg Dref Iref IrefGray
            
            [Row Col K] = size(S1);
            N           = Row*Col;
            
            %%Shape B
            shapeB  = gampdf(1:K,2,30)* K; 
            shapeB  = repmat(shapeB(:)',N,1);
            shapeB  = reshape(shapeB,[Row,Col,K]);
            
            %% True IRF
            load F_real2_100s
            F = F(:,100); F(F<(0.01*max(F(:))))=0;
            h = [zeros(100,1); F(2:2:end); zeros(100,1)];
            h  = h/sum(h);
            [vvv, MaxF] = max(h); clear F
            for i=1:length(h)
                F(:,i)  = circshift(h,[i-MaxF 0]);
            end
            F = F(1:K,1:K);
            
            %% SBR, PPP levels
%             pppY = logspace(log10(0.1),log10(1000),ppplength);
%             SBR  = logspace(log10(0.01),log10(100),ppplength);
            pppY = [1.0 4.0 16.0 64.0];
            SBR  = [1.0 4.0 16.0 64.0];
            
            %% Convolution in Times by IRF
            S1_Conv =  reshape((F*reshape(S1(:,:,:),[N K])')', [Row Col K]) ;
            S2_Conv =  reshape((F*reshape(S2(:,:,:),[N K])')', [Row Col K]) ;
            S3_Conv =  reshape((F*reshape(S3(:,:,:),[N K])')', [Row Col K]) ;
            Sg_Conv =  reshape((F*reshape(Sg(:,:,:),[N K])')', [Row Col K]) ;
            
            for LevP = 1:length(pppY)
                LevP
                for LevSBR = 1:length(SBR)
                    Lev_S = SBR(LevSBR) .*pppY(LevP) ./(1+SBR(LevSBR));
                    Lev_B = pppY(LevP)- Lev_S;
                    SBRmat(LevP,LevSBR)  = SBR(LevSBR);
                    pppYmat(LevP,LevSBR) = pppY(LevP);
                    %                     S1_GT(:,:,:,LevP,LevSBR)   = S1*Lev_S;
                    %                     S2_GT(:,:,:,LevP,LevSBR)   = S2*Lev_S;
                    %                     S3_GT(:,:,:,LevP,LevSBR)   = S3*Lev_S;
                    %                     Sg_GT(:,:,:,LevP,LevSBR)   = Sg*Lev_S;
                    S1_Poiss(:,:,:,LevP,LevSBR) = poissrnd(S1_Conv*Lev_S);
                    S2_Poiss(:,:,:,LevP,LevSBR) = poissrnd(S2_Conv*Lev_S);
                    S3_Poiss(:,:,:,LevP,LevSBR) = poissrnd(S3_Conv*Lev_S);
                    Sg_Poiss(:,:,:,LevP,LevSBR) = poissrnd(Sg_Conv*Lev_S);
                    
                    Y1_Poiss(:,:,:,LevP,LevSBR) =  S1_Poiss(:,:,:,LevP,LevSBR) + poissrnd(Lev_B*shapeB/K);
                    Y2_Poiss(:,:,:,LevP,LevSBR) =  S2_Poiss(:,:,:,LevP,LevSBR) + poissrnd(Lev_B*shapeB/K);
                    Y3_Poiss(:,:,:,LevP,LevSBR) =  S3_Poiss(:,:,:,LevP,LevSBR) + poissrnd(Lev_B*shapeB/K);
%                     Yg_Poiss(:,:,:,LevP,LevSBR) =  Sg_Poiss(:,:,:,LevP,LevSBR) + poissrnd(Lev_B*shapeB/K);
                    Yg =  Sg_Poiss(:,:,:,LevP,LevSBR) + poissrnd(Lev_B*shapeB/K);
        
                    save(['/Users/jk2021/work/unroll_plot/data/middlebury1024gamma/' + string(Lev_P) + string(Lev_SBR) + '.mat'], 'Y')
                end
            end
%             Name = ['Art_Noisy_3Wav_1peak_185_232_300_5PPPy_5SBR_GammaBack.mat'];
        

%             Name = ['Art_Noisy_3Wav_1peak_185_232_300_5PPPy_5SBR_GammaBack_IndD.mat'];
            
%             save(Name,'SBRmat','pppYmat','pppY','SBR', 'F',  'Dref', 'Iref', 'IrefGray',...
%                 'S1','S2','S3','Sg', ...
%                 'Y1_Poiss','Y2_Poiss','Y3_Poiss','Yg_Poiss',...
%                 'S1_Poiss','S2_Poiss','S3_Poiss','Sg_Poiss','-v7.3')
             
    end
end
