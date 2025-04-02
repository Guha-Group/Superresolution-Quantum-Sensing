function q0 = SimulateBSPADE(x0,s,kappa,n,sigma)
    % point source parameters
    x1 = x0-s; b1 = .5-kappa;
    x2 = x0+s; b2 = .5+kappa;

    % probability of detection in HG00 mode
    p0 = b1*poisspdf(0,(x1/2*sigma)^2) + b2*poisspdf(0,(x2/2*sigma)^2);
    
    % sample total fraction of photons in HG00 mode from binomial
    q0 = binornd(n,p0,1);
end