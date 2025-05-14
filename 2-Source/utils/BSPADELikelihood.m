function L_BS = BSPADELikelihood(q0,n,x0,s,kappa,sigma,norm_approx)
    % Computes the binomial likelihood for binary SPADE in the case where
    % q0 photons were observed in the fundamental mode. The probability of
    % photon detection in the fundamental mode p0 is function of
    % x0,s,kappa. We use the normal approximation of the binomial wherever
    % applicable to speed up computations for large array queries.
    %
    % q0:   number of photons observed in fundamental mode
    % n:    number of photons detected in total
    % x0:   midpoint of point sources
    % s:    half-separation of point sources
    % kappa:brightness bias between point sources

    % point source parameters
    x1 = x0-s; b1 = .5-kappa;
    x2 = x0+s; b2 = .5+kappa;
    
    % probability of detection in HG00 mode
    p0 = b1*exp(-(x1/2*sigma).^2) + b2*exp(-(x2/2*sigma).^2);

    % promote array dimensions
    Q0 = repmat(q0,size(p0)); P0 = repmat(p0,size(q0));

    % compute the BSPADE likelihood
    if norm_approx
        %criterion for normal approximation to binomial
        use_normal_approx = n.*p0.*(1-p0)>20;
        L_BS = normpdf(q0,n*p0,sqrt(n.*p0.*(1-p0)));
        L_BS(:,~use_normal_approx) = binopdf(Q0(:,~use_normal_approx),n,P0(:,~use_normal_approx));
    else
        L_BS = binopdf(Q0,n,P0);
    end
end


