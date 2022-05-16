function [Y2] = Duplicate_Borders(Y,DupPix);
 
[r c T]   = size(Y);
Y2        = zeros(r+2*DupPix,c+2*DupPix,T);    
[r2 c2 T] = size(Y2);


Y2(DupPix+1:r2-DupPix,DupPix+1:c2-DupPix,:)  = Y;

% Y2(1:DupPix,DupPix+1:c2-DupPix,:)            = Y(1:DupPix,:,:);
% Y2(r+DupPix:r+2*DupPix,DupPix+1:c2-DupPix,:)  = Y(r:-1:r-DupPix,:,:);
% 
% Y2(:,DupPix:-1:1,:)  = Y2(:,DupPix+1:DupPix+DupPix,:);
% Y2(:,c+DupPix:c+2*DupPix,:)  = Y2(:,c+DupPix:-1:c,:);

  
Y2(1:DupPix,DupPix+1:c2-DupPix,:)               = Y2(2*DupPix:-1:DupPix+1,DupPix+1:c2-DupPix,:);
Y2(r+DupPix+1:r+2*DupPix,DupPix+1:c2-DupPix,:)  = Y2(r+DupPix:-1:r+1,DupPix+1:c2-DupPix,:) ;
Y2(:,DupPix:-1:1,:)            = Y2(:,DupPix+1:2*DupPix,:);
Y2(:,c+DupPix+1:c+2*DupPix,:)  = Y2(:,c+DupPix:-1:c+1,:);

 

