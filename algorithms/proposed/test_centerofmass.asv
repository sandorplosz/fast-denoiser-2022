T=10;
N=1;
window=[1,1];

Dep=zeros(N, 1);
Yt=randi(50,T,N);
Yt=Yt-20;
Yt(Yt<0)=0;
for i=1:T
    fprintf("%i, ", Yt(i));
end
xmax = randi(T,1);   

for i=1:N
    v=max(xmax(i)-window(1),1):min(xmax(i)+window(2),T);
    Dep(i)=v*Yt(v,i)/sum(eps+Yt(v,i))
end


%%%%%%%%%%%%%%%
for xmax=1:T 
    i=1;
    v=max(xmax(i)-window(1),1):min(xmax(i)+window(2),T);
    Dep(i)=v*Yt(v,i)/sum(eps+Yt(v,i));
    fprintf("%f, ", Dep);
end

