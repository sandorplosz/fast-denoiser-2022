function out = fcn_convolve(Y,h,w)

T = length(Y);


d = ifft(Y.*fft(log(w.*h*T+1)),'symmetric');

d = d(round(T/4):round(T*3/4),:);

offset = max(d);
out = log(mean(exp(d-offset)))+offset;

%disp(length(w));

% d2 = log(T*w+1).*ifft(Y.*fft(h),'symmetric');
% out2 = log(mean(exp(d2-offset)))+offset;

end