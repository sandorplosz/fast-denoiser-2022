function [det,w_map] = detect(y,h,max_r)


if sum(y)==0
    det = 0;
    w_map = 0;
else
    T = length(y);
    
    %% priors
    %SBR = 0.5;
    alpha_r = 2;
    beta_r = alpha_r/(max_r/2);
    
    sumy = sum(y);
    
    alpha_b = 1;
    beta_b = alpha_b/(max_r/T);
    
    %c0 = alpha_r*log(T*beta_r)+alpha_b*log(beta_b)-gammaln(alpha_r)-gammaln(alpha_b)+gammaln(alpha_b+alpha_r);
    %prior_w  =@(w) (alpha_r-1)*log(w)-(alpha_r+alpha_b)*log(beta_b+beta_r*T*w)+c0;
    
    h = flipud(h);
    Y = fft(y);
    
    %% densities
    c1 = gammaln(sumy+alpha_b+alpha_r)-gammaln(alpha_r)+alpha_r*log(T*beta_r);
    p0 = -(sumy+alpha_b)*log(T+beta_b)+gammaln(sumy+alpha_b);
    
    p1 = @(w) (alpha_r-1)*log(w)+fcn_convolve(Y,h,w)-(sumy+alpha_r+alpha_b)*log(T*(1+w*(1+beta_r))+beta_b)+ c1;
    
    %% find map
    w_int = logspace(-2,2,20);
    out = p1(w_int);
    [offset,ind]=max(out);
    w_map = w_int(ind);
   
    
    %% output 
    dens_cero = p0-offset;
    fun = @(w) exp(p1(w)-offset);
    tol = exp(dens_cero)/10;
    integ = log(integral(fun,0,Inf,'AbsTol',tol));
    det = integ-dens_cero;
    
    if isnan(det)
        keyboard
    end
end

end