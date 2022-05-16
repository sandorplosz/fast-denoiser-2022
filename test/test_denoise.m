clearvars
pars.sigIRF=5;

dn_orig = randi (999, 1, 25);
dn_orig = max(1,(dn_orig-300)/10);
rn_orig = randi(60, 1, 25);
rn_orig = max(0, rn_orig-20);

%dn = [0.0 38.0 35.5 0.0 0.0 19.8 65.9 4.1 28.5 0.0 45.1 0.0 20.6 39.9 59.1 65.9 24.7 0.0 0.0 0.0 54.0 0.0 51.4 0.0 62.9];
%rn = [0 28 0 12 0 17 0 20 22 25 8 0 0 35 0 30 13 40 0 7 0 38 0 27 30];

for i = 1:length(dn_orig)
    fprintf("%.1ff, ", dn_orig(i));
    if mod(i,10)==0
        fprintf("\n");
    end
end
fprintf("\n");
for i = 1:length(rn_orig)
    fprintf("%i, ", rn_orig(i));
end
fprintf("\n");

dn=dn_orig;
rn=rn_orig;

a=min(mean(rn(dn>1)), mean(rn(dn==1)));
rn(dn==1) = a;
nnd = 1; dim=1;

mask = rn==0;
dn(mask) = [];
rn(mask) = [];        
if(isempty(dn))            
    DhrNew3(nnd) = 1;            
    IntTemp(nnd,dim) =0;        
else           
    DiffV = abs(dn - dn(:));            
    vvvb =  (rn(:)*rn) .* exp( - DiffV.^2/(2*pars.sigIRF^2));            
    [~, indD] = max(sum(vvvb));     %********            
    %DhrNew3(nnd) = dn(indD);     %********            
    Int = sum((DiffV(:,indD)<(pars.sigIRF*3)) .*rn(:));          
    DhrNew3 = sum((DiffV(:,indD)<(pars.sigIRF*3)) .*dn(:).*rn(:))/ (eps+sum((DiffV(:,indD)<(pars.sigIRF*3)) .*rn(:)));        
end  

Int = Int / length(dn_orig);
