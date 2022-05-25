downSam = 2;
K = 1024;
s_back = ["UnifBack", "GammaBack"];
s_backest=["nobackest", "withbackest1", "withbackest2", "withbackest3"];    

algodir = mfilename('fullpath');
algodir = erase(algodir, "init");
addpath([algodir, '/proposed/']);
basedir = [algodir, '/../'];
%dir(basedir)

scenes = { 'Reindeer', 'Art', 'Plastic', 'Moebius', 'Laundry', ...
            'Dolls', 'Bowling1', 'Books' };
selectedScene=scenes{2};

irfs = ["F_gauss_sig_6", "F_real_proc"];

load(strcat(basedir, 'data_generation/',irfs(selirf),'.mat'))
h=F(:,ceil(size(F,2)/2));
F_window=getIRFWindow(h,0.1);
%F = processF(F,K);
load([basedir, 'data_generation/Art_ref_img.mat']);

dataDir = strcat(basedir, 'data_generation/', irfs(selirf));
[row, col] = size(Dref);

params = loadParams(row*col, F_window);
Dref_orig = Dref;
Dref = Dref*params.Tbin*3*10^8/2;