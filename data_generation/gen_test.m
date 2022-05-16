T=100;
PPP=20;
sigma2=3^2;
h= 1/sqrt(2*pi*sigma2)* exp(- ((1:T)-round(T/2)).^2/2/(sigma2));
h = h(:)/sum(h);
[~,attack] = max(h);
l = length(h);
F=zeros(l,l);    
for i=1:l
    F(:,i)  = circshift(h,[i-attack 0]);
end

a=randi(T,1,1);
b=zeros(1,T);
b(a)=PPP;
s=F*b';
c=poissrnd(s);
sum(c)
%figure; plot(c)

Tof=[];
ind = find(c>0);
for k=1:length(ind)
    Tof = [Tof kron(ind(k),ones(1,c(ind(k))) )];
end
%Tof    = Tof(randperm(length(Tof)));  

vs=length(Tof);
Tofp = [Tof zeros(1,32-vs)];
for i=1:length(Tofp)
    fprintf("%i, ", Tofp(i));
end

detect_AH_NonUnifBack_v7_v2(Tof ,sigma2, 1, 1, T, 10 , 0.5, 252, 3, 1/T*ones(1,T), 0);