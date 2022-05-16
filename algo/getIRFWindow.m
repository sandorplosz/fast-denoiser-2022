function window = getIRFWindow(h, thresh)
    [m,ind]=max(h);
    mid=ceil(length(h)/2);
    h=circshift(h,[mid-ind 0]);
    window(1)=sum(h(1:mid)>thresh*m);
    window(2)=sum(h(mid:end)>thresh*m);
end