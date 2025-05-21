function x = dirchrnd(alpha)
    % randomly samples dirichlet random variables with rate parameters
    % given by alpha [K,1]
    y = gamrnd(1,alpha);
    x = y./sum(y,1);
end