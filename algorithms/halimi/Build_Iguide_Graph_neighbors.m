function [ I_up_bar, Iguide, LvectSR ]   = Build_Iguide_Graph_neighbors(InMat,r,c,NeighboursSR,Neighbours,GuideI);


LvectSR      = size(NeighboursSR,2);
[N,L,R0]     = size(InMat);

I_up_bar = [];pp=1;
for ell =1:R0
    for k=1:LvectSR
        for wavelength = 1:L
            I_up_bar(pp,:,wavelength) =   InMat(NeighboursSR(:,k),wavelength,ell)  ;% nd SCale  N wavelength
        end
        pp = pp+1;
    end
end

switch GuideI
    case 0
        Iguide = I_up_bar;%.*(Dguide>0);
        
    case 1
        
        II  = permute(I_up_bar,[2,1,3]); % N 27 3
        % % % %%% No filtering
%         keyboard
 
        
        ChooseMeth =  1%(mean(mean(II(:,1:LvectSR:end,:)))>3);
        if(ChooseMeth)
            im   = reshape(II(:,1:LvectSR:end,:),[r,c,R0,L])   ; % r c Scales  wavelength
            IntR = 1:R0;
        else
            im  = mean(reshape(II(:,1:LvectSR:end,:),[r,c,R0,L]),3)   ; % r c Scales  wavelength %%%%% Method1
            IntR = 1;
        end
%         %%%%%% DnCNN
%         net = denoisingNetwork('DnCNN');
%         for ell = IntR       
%             for Wave=1:L
%                 FiltI(:,:,ell,Wave) = denoiseImage(im(:,:,ell,Wave),net);%denoiseImage(Ymat2(:,:,r,1),net);%
%             end
%         end
        %%%%%% VST
        cd E:\Travail_30_09_2015\ahalimi\Post_doc1\Codes\iterVSTpoisson_STANDALONE    
        for ell = IntR      %%%%% Method1
            for Wave=1:L
                FiltI(:,:,ell,Wave)=iterVSTpoisson(im(:,:,ell,Wave));            
            end
        end
        cd E:\Travail_30_09_2015\ahalimi\Post_doc1\PC_Lidar_Graph\Main_FastDenoiseFog_1peak\Reconstruction_Bayesian_L1\Demo_Denoising
%         keyboard
        if(ChooseMeth==0)
            FiltI = repmat(FiltI,[1, 1,R0,1]); %%%%% Method1
        end
        FiltI = reshape(FiltI, [N, R0, L]) ; %  N R L
        Iguide = [];pp=1;
        for ell =1:R0
            for k=1:LvectSR
                for wavelength = 1:L
                    Iguide(pp,:,wavelength) =   FiltI(NeighboursSR(:,k),ell,wavelength)  ;% nd SCale  N wavelength
                end
                pp = pp+1;
            end
        end
end




