function p = loadParams(N, F_window)

    p = struct(...
        'Alpha',        1,  ...
        'Beta',         1,  ...  % non-informative prior 
        'NbrePhoton',   10, ...  % Approximation with 10 photons (choose <12); 
        'ProbPrior',    0.5*ones(N,1), ...
        'Tbin',         20*10^(-12), ... % time sample or bin in seconds
        'SLight',       3*10^8, ...        
        'CvPC',         0 ...      % 1:conv,   0: operations on PC    
    );

    p.limitC = nchoosek(p.NbrePhoton,round(p.NbrePhoton/2));
    p.Attack = F_window(1);
    p.trailing = F_window(2);
    p.sigIRF=round(min(F_window)*2/4.3,1);
    p.IRFw = p.Attack+p.trailing-1;
    p.ThreshDep = (-9/198*p.Tbin*10^(12) + 10+18/198)*p.SLight/2*p.Tbin;

end