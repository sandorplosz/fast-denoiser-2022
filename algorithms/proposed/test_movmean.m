clearvars;
a=randi(20,1,100);
a=max(a-3,0);
w=[3,2];

ref=movmean(a,w);

% Solution 1
% We pad the array with 0 to be exclusive sum
b=[0, cumsum(a)];
%temp=[[0, a];b];
l=length(b);
r=zeros(1,l);
for i=1:l
    w1=max(i-w(1)-1,1);
    w2=min(i+w(2),l);    
    r(i)=(b(w2)-b(w1))/(w2-w1);
end
r=r(2:l);
isequal(r,ref)


% Solution 2
% No padding, but additional check
b=cumsum(a);
l=length(b);
temp=[a;b];
r=zeros(1,l);
for i=1:l
    w1=max(i-w(1)-1,0);
    w2=min(i+w(2),l);    
    if(w1>0)
        v=b(w1);
    else
        v=0;
    end
    r(i)=(b(w2)-v)/(w2-w1);
end
isequal(r,ref)