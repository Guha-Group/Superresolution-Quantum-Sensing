function x = SimulateDirectImaging(x0,s,kappa,n,sigma)
    % point source parameters
    x1 = x0-s; b1 = .5-kappa;
    x2 = x0+s; %b2 = .5+kappa;
    
    n1 = binornd(n,b1,1); % photons coming from first source
    n2 = n-n1;            % photons coming from second source
    x = [normrnd(x1,sigma,n1,1); normrnd(x2,sigma,[n2,1])];
end

