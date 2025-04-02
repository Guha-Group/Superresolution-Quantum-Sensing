function L_BS = BSPADELikelihood(q0,n,x0,s,kappa,sigma)
    % point source parameters
    x1 = x0-s; b1 = .5-kappa;
    x2 = x0+s; b2 = .5+kappa;
    
    % probability of detection in HG00 mode
%   p0 = b1*poisspdf(0,(x1/2*sigma).^2) + b2*poisspdf(0,(x2/2*sigma).^2);
    p0 = b1*exp(-(x1/2*sigma).^2) + b2*exp(-(x2/2*sigma).^2);
    L_BS = binopdf(repmat(q0,size(p0)),n,repmat(p0,size(q0)));
end