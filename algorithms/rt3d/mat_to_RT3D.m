function  mat_to_RT3D(Y, h, scale_ratio, filename)

%% FILE STRUCTURE
% Nrow(uint16) Ncol(uint16) T(uint16) scale_ratio(float) ManyIRFflag(uint16)
%impulse_len(uint16)  h1(float)...hN(float)  bin_index{1}bin_counts{1}0xFFFF bin_index{2}bin_counts{2}0xFFFF ... 0xFFFF bin_index{NrNc}bin_counts{NrNc} 
ordering = 'n';

%filename = erase(filename,".mat");

%% load dataset
[Nrow,Ncol,T,frames] = size(Y);


%% crop impulse response (remove very low values)
if numel(h)>size(h,1)
    min_ind_p = T;
    pixel_h = [1,1];
    for j=1:Ncol
        for i=1:Nrow
            [~,ind_p] = max(h(i,j,:));
            if ind_p<min_ind_p
                min_ind_p = ind_p;
                pixel_h = [i,j];
            end
        end
    end
    r = squeeze(h(pixel_h(1),pixel_h(2),:));
else
    r = h;
end

r_max = max(r);
thres = r_max*0.01;
t = 1;
while r(t)<thres
    t = t + 1;
end
start_t = max([t-1,1]);
t = T;
while r(t)<thres 
    t = t - 1;
end
end_t = t;

if numel(h)>size(h,1)
    h = h(:,:,start_t:end_t);
else
    h = h(start_t:end_t);
end

%% prepare RT3D dataset file
filename = [filename '.rbin'];
    
fileID = fopen(filename,'w');

fwrite(fileID,Nrow,'uint16',ordering);
fwrite(fileID,Ncol,'uint16',ordering);
fwrite(fileID,T,'uint16',ordering);
fwrite(fileID,scale_ratio,'float',ordering);

%% save impulse response
if numel(h)>size(h,1)
    % multiple IRFs (one per pixel)
    h = round(h/max(h(:))*(2^16-1));
    fwrite(fileID,hex2dec('FFFF'),'uint16',ordering);
    fwrite(fileID,size(h,3),'uint16',ordering);
    for j=1:Ncol
        for i=1:Nrow
            for t=1:size(h,3)
                fwrite(fileID,h(i,j,t),'uint16',ordering);
            end
        end
    end
else
    % single IRF
    fwrite(fileID,hex2dec('0000'),'uint16',ordering);
    fwrite(fileID,length(h),'uint16',ordering);
    for i=1:length(h)
        fwrite(fileID,h(i),'float',ordering);
    end
end


%% save photon detections
fwrite(fileID,frames,'uint16',ordering);

for f=1:frames
    disp(['Converting frame ' num2str(f) ' out of ' num2str(frames)])
    Z = reshape(Y(:,:,:,f),Nrow*Ncol,T);
    for n=1:size(Z,1)
        
        z = find(Z(n,:)>0);
        for i=1:length(z)
            fwrite(fileID,z(i),'uint16',ordering);
            fwrite(fileID,Z(n,z(i)),'uint16',ordering); 
        end
        fwrite(fileID,hex2dec('FFFF'),'uint16',ordering);
        
    end
end

%% close dataset file
fclose(fileID);
disp('Done!')

end

