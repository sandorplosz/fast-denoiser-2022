%% Generate synthetic data - TCI Rev1
%% this is the data used for the TCI paper after revision
%% by Abderrahim Halimi
clear all
clc
choice = 1; %0  one peak,   1: multipeaks
D_HR   = double(imread('disp1.png'));
I_HR   = double(imread('view1.png'));
downSam = 2;

Choice{1}  =   'GENERATE_CLEAN_CUBES';%   
Choice{2}  =   'GENERATE_NOISY_CUBES';%      Unif Back
Choice{3}  =   'GENERATE_NOISY_CUBES_GAMMA';%     Gamma back 

for Indchoix = 2
    choix = Choice{Indchoix};
    switch choix
        case 'GENERATE_CLEAN_CUBES'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Median filter: Reconstruct missing parts
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [Row Col]    = size(D_HR);
            K            = 300;
            [indI  indJ] = find(D_HR==0);
            DiffD = (min(cat(3,abs(D_HR - circshift(D_HR,[-1,0])),abs(D_HR - circshift(D_HR,[1,0]))...
                ,abs(D_HR - circshift(D_HR,[0,1])) , abs(D_HR - circshift(D_HR,[0,-1]))),[],3));
            [indI1  indJ1] = find(DiffD > 50);
            [indI  ] = [indI ; indI1];
            [indJ  ] = [indJ ; indJ1];
            
            step         = 10;
            D0_HR        = D_HR;
            for n=1:length(indI)
                xi  = max(1,indI(n)-step):min(Row,indI(n)+step);
                yj  = max(1,indJ(n)-step):min(Col,indJ(n)+step);
                vvv = D_HR(xi,yj);
                D0_HR(indI(n),indJ(n)) = round(median(vvv(:)));
            end
            figure
            subplot(2,1,1);imagesc(D_HR)
            subplot(2,1,2);imagesc(D0_HR)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Downsampling
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Dref     =  max(1,D0_HR(2:downSam:end,1:downSam:end));
            Iref     =  max(1,I_HR(2:downSam:end,1:downSam:end,:));
            
            [Row Col]  = size(Dref);
            %            Iref = Iref/sum(Iref(:));
            Iref = Iref/max(Iref(:));
            
            IrefGray = rgb2gray(Iref);  
            IrefGray = IrefGray/mean(IrefGray(:));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% 3D data cube: 3 wavelengths
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             S1 = zeros(Row,Col,K);
%             S2 = zeros(Row,Col,K);
%             S3 = zeros(Row,Col,K);
            Sg = zeros(Row,Col,K);
            Dref = K-Dref;

            for i=1:Row
                i
                for j=1:Col
%                     S1(i,j,Dref(i,j)) = (Iref(i,j,1));
%                     S2(i,j,Dref(i,j)) = (Iref(i,j,2));
%                     S3(i,j,Dref(i,j)) = (Iref(i,j,3));
                    Sg(i,j,Dref(i,j)) = (IrefGray(i,j));
                end
            end
%             figure;imagesc(sum(S1,3))
            
            figure; I1=IrefGray(:,:,1);
            scatter3(kron([1:Col]',ones(Row,1)),Dref(:),kron(ones(Col,1),[Row:-1:1]'),1,I1(:),'.'),
            [Row Col K] = size(Sg);
            
%             figure;imagesc(sum(S1~=0,3));
            %
            %  figure; I1=I(:,:,1);
            %  scatter3(kron([1:Col]',ones(Row,1)),Dref(:),kron(ones(Col,1),[Row:-1:1]'),1,I1(:),'.'),
            
%             save Art_GT_1Wavelengths_1peak_185_232_300.mat S1 S2 S3 Sg Dref Iref IrefGray
            save Art_GT_1Wavelengths_1peak_185_232_300_IndD.mat Sg Dref Iref IrefGray
           
        case 'GENERATE_NOISY_CUBES'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%  Noise
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% Generate noisy data/different Counts per pixel
%             load Art_GT_1Wavelengths_1peak_185_232_300.mat S1 S2 S3 Sg Dref Iref IrefGray
            load Art_GT_1Wavelengths_1peak_185_232_300_IndD.mat Sg Dref Iref IrefGray

            [Row Col, K] = size(Sg);
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
            pppY = logspace(log10(0.1),log10(1000),5);
            SBR  = logspace(log10(0.01),log10(100),5);
            
            %% Convolution in Times by IRF
%             S1_Conv =  reshape((F*reshape(S1(:,:,:),[N K])')', [Row Col K]) ;
%             S2_Conv =  reshape((F*reshape(S2(:,:,:),[N K])')', [Row Col K]) ;
%             S3_Conv =  reshape((F*reshape(S3(:,:,:),[N K])')', [Row Col K]) ;
            Sg_Conv =  reshape((F*reshape(Sg(:,:,:),[N K])')', [Row Col K]) ;
            
            for LevP = 1:length(pppY)
                LevP
                for LevSBR = 1:length(SBR)
                    Lev_S = SBR(LevSBR) .*pppY(LevP) ./(1+SBR(LevSBR));
                    Lev_B = pppY(LevP)- Lev_S;
                    SBRmat(LevP,LevSBR)  = SBR(LevSBR);
                    pppYmat(LevP,LevSBR) = pppY(LevP);
                    
                    Sg_Poiss(:,:,:,LevP,LevSBR) = poissrnd(Sg_Conv*Lev_S);
                    
                    Yg_Poiss(:,:,:,LevP,LevSBR) =  Sg_Poiss(:,:,:,LevP,LevSBR) + poissrnd(Lev_B*ones(Row, Col, K)/K);
                    
                end
            end
              Name = ['Art_Noisy_1Wav_1peak_185_232_300_5PPPy_5SBR_UnifBack_IndD.mat'];

            save(Name,'SBRmat','pppYmat','pppY','SBR', 'F',  'Dref', 'Iref', 'IrefGray',...
                'S1','S2','S3','Sg', ...
                'Yg_Poiss',...
                'Sg_Poiss','-v7.3')
    end
end
