%% Constants
L=1;
SLight     =  3*10^8 ;
FramHist   = 0;

%% Data or system related constants 
ind3D1   = kron((1:col)', ones(row,1));
ind3D3   = kron(ones(col,1),(row:-1:1)');
WindSR   = 3; % Size of window for ...
WindMed  = 5; % Size of window for median
Attack   = F_window(1) *ones(1,L);%7  ; % left side of IRF
trailing = F_window(2)*ones(1,L);%50;  % right side of IRF
IRFw        = 2*min([Attack,trailing]); % width of IRF
I_resol     = [1 3 7] ;  % size of spatial correlations Requires 1x1  3x3 9x9
Tbin        = 20*10^(-12); % time sample or bin in seconds
ThreshDep   = (-9/198*Tbin*10^(12) + 10+18/198)  *SLight/2*Tbin; %1bin for 200ps, 10bins for 2ps
sigIRF      = round(min(F_window)*2/4.3,1);
WindWM      = 9;   % Size of window for weighted median
weightCurr  = 0.3; % Weight for current sample compared to past
Dref        = Dref  *SLight/2*Tbin;

%% Computed constants / Graph of pixel correlations
N             = row*col;
[Neighbours_Med indGraph_Med,local_shifts_Med]  =  Build_Graph_Neighbours_Array_v2(row,col,[1 WindMed]); % Define graph of correlations between pixels
[Neighbours_WM, indGraph_WM,local_shifts_WM]    =  Build_Graph_Neighbours_Array_v2(row,col,[1 WindWM]);  % Define graph of correlations between pixels
Neighbours_Med = Neighbours_Med(1:end-1,:);
Neighbours_WM  = Neighbours_WM(1:end-1,:);
nd             = size(local_shifts_WM,1);
weights        = [1./sum(abs(local_shifts_WM(2:end,:).^2)') ];
for t=1:FramHist
    weights = [weights 1./(t+sum(abs(local_shifts_WM.^2)'))];
end
weights        = [weightCurr weights/sum(weights(:))*(1-weightCurr)]; % 1  L-1  T{L L L}  x  N = TL x N

SR =0;% 1;
%     keyboard
rowSR = (row*(2^SR )-(2^SR -1));
colSR = (col*(2^SR )-(2^SR -1));
[Neighbours_MedSR indGraph_MedSR,local_shifts_MedSR] =  Build_Graph_Neighbours_Array_v2(rowSR,colSR,[1 floor(WindSR/2)*((2^SR -1)+1)*2+1]);% Define graph of correlations between pixels
Neighbours_MedSR  = Neighbours_MedSR(1:end-1,:);

mask        = ones(row,col); 