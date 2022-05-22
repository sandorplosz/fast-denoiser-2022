useTargetDetect = 0; 
% F_gauss_sig_6=1,  F_real_proc=2
selirf=2;
run ../algorithms/proposed/init.m

myFolder = strcat('./', irfs(selirf), '/');
myFiles = dir(strcat(myFolder, 'Samples*.mat'));

maxP=0;

for k = 1:length(myFiles)
  fname = myFiles(k).name;
  load([myFiles(k).folder, '/', fname]);
  pv = sum(Y,2);
  mx = full(max(pv));
  fprintf(1, 'File %s maxPh=%i, meanPh=%f\n', fname, mx, full(mean(pv)));
  if(mx>maxP), maxP=max(pv); end
end