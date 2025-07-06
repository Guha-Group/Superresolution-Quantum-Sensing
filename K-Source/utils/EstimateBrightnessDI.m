function b_est_prev = EstimateBrightnessDI(xy_samples,xy,sigma)
    % Estimates the brightnesses of various point sources from a direct 
    % imaging measurement. The positions of the sources are assumed to be
    % known and given by xy.

    % initialize the estimate
    num_sources = size(xy,1);
    b_est_prev = ones([num_sources,1])/num_sources;

    % conditional probability of photon arrival location given emitter
    % position mu.
    W = zeros(size(xy_samples,1),num_sources);
    for k = 1:num_sources
        W(:,k) = mvnpdf(xy_samples,xy(k,:),eye(2)*sigma^2);
    end

    % iteratively update the brightness estimates (maximum likelihood)
    converged = 0;
    max_iters = 1e3;
    t=0;
    while ~converged

        % compute soft assignment weights
        weights = W .* b_est_prev';
        weights = weights ./ sum(weights,2);

        % update the brightness estimate
        b_est = mean(weights,1)';
        
        % check convergence criteria
        converged = (norm(b_est-b_est_prev)<1e-3) || (t==max_iters);

        % update previous estimate
        b_est_prev = b_est;

        % update iteration counter
        t = t+1;
    end
end