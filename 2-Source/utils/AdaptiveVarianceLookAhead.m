function E_var_s = AdaptiveVarianceLookAhead(M,m1,ps,pe,s,e,sigma)
    m2 = M-m1;
    q = (0:m2)';
    ds = s(2)-s(1);
    de = e(2)-e(1);

    % compute the look-ahead variances for each possible value of BSPADE outcome 
    L_BS = BSPADELikelihood(q,m2,e,s,0,sigma);          % BSPADE likelihood
    p_sIeq = (L_BS.*ps)./sum(L_BS.*ps*ds,3);            % posterior probability on separation
    p_sIq = sum(p_sIeq.*pe*de,2);                       % posterior probability on separation marginalized over the pointing error
    var_sIq = sum(s.^2.*p_sIq*ds,3) - sum(s.*p_sIq*ds,3).^2;    % variance of separation under the posterior for all possible BSPADE measurement outcomes

    % compute the expectation of the look-ahead variances under the current
    % separation estimate
    s_hat = sum(s.*ps*ds,2);
    weights = sum(BSPADELikelihood(q,m2,e,s_hat,0,sigma).*pe*de,2);
    E_var_s = sum(weights.*var_sIq,1);
end


