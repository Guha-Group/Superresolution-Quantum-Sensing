function L_DI_alt = DirectImagingLikelihood_approx(x,x0,s,kappa,sigma,prop_flag)
    % computes an approximate likelihood for direct imaging under the
    % approximation that the source separation is small. If prop_flag is
    % on, just return a value proportional to the likelihood (which is more
    % numerically stable) that ignores the common PSF which is independent
    % of kappa.
    
    L_DI_alt = ones([1,numel(x0),numel(s),numel(kappa)]);
    for j = 1:numel(x)
        L_DI_alt = L_DI_alt .* (1+2*kappa.*s.*(x(j)-x0)/(sigma^2));
    end
    
    if ~prop_flag
        % include common gaussian factors
        gaussfactor = prod(normpdf(x,x0,sigma),1);
        L_DI_alt = gaussfactor*L_DI_alt;
    end 
end