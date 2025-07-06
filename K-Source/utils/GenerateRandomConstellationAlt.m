function xyb = GenerateRandomConstellationAlt(num_sources,min_sep_frac,alpha_vec,sigma)
    % generates a random constellation where all the sources are separated
    % by at least min_sep_frac. The brightnesses are sampled from a
    % Dirichlet distribution with parameters alpha.
    xyb = GenerateRandomConstellation(num_sources,min_sep_frac,1,sigma);
    b = dirchrnd(alpha_vec);
    xyb(:,3) = b;
end