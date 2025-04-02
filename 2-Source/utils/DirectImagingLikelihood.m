function L_DI = DirectImagingLikelihood(x,x0,s,kappa,sigma)
    % point source parameters
    x1 = x0-s; b1 = .5-kappa;
    x2 = x0+s; b2 = .5+kappa;
    
    L_DI = ones([1,numel(x0),numel(s),numel(kappa)]);
    %logL_DI = zeros(size(L_DI));
    for j = 1:numel(x)
        L_DI = L_DI .* (b1.*normpdf(x(j),x1,sigma)+b2.*normpdf(x(j),x2,sigma));
        %logL_DI = logL_DI + log(b1.*normpdf(x(j),x1,sigma)+b2.*normpdf(x(j),x2,sigma));
    end
end