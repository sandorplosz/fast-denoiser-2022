for i=1:200
    n=9;
    if i==1
        dn=zeros(1,n);
    else    
        dn=randi(5,1,n)-1;
    end
    rn=randi(5,1,n)-1;    
    rn(dn==1)=0;
    mask = dn==0;      
    
    % Reference
    dn1 = dn(~mask);
    rn1 = rn(~mask);  
    if(isempty(dn1))            
        DhrNew1 = 1;            
        IntNew1 =0;      
    else
        DiffV = abs(dn1 - dn1(:));            
        vvvb1 =  (rn1(:)*rn1) .* exp( - DiffV.^2/(2*pars.sigIRF^2));     
        [~, indD] = max(sum(vvvb1));     %********            
        IntNew1 = sum((DiffV(:,indD)<(pars.sigIRF*3)) .*rn1(:));
        DhrNew1 = sum((DiffV(:,indD)<(pars.sigIRF*3)) .*dn1(:).*rn1(:))/ (eps+sum((DiffV(:,indD)<(pars.sigIRF*3)) .*rn1(:)));
    end
    
    % CUDA adapted version
    DiffV2 = abs(dn - dn(:));            
    vvvb2 =  (rn(:)*rn) .* exp( - DiffV2.^2/(2*pars.sigIRF^2));   
    vvvb2(:,mask)=0;
    vvvb2(mask,:)=0;    
    [~, indD] = max(sum(vvvb2));     %********     
    mask2 = ~mask(:) & DiffV2(:,indD)<(pars.sigIRF*3);
    IntNew2 = sum(mask2.*rn(:));
    DhrNew2 = sum(mask2.*dn(:).*rn(:))/ (eps+sum(mask2.*rn(:))); 
    if DhrNew2==0
        DhrNew2=1;
    end
    
    if(DhrNew1~=DhrNew2)
        fprintf("NEQ\n");
        break;
    end
end
