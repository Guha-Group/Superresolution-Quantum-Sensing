function xy = DILocalizeSources(xy0,xy_samples,sigma)
% Localize sources from direct imaging measurement. Uses expectaction
% maximization under the assumption that the source brightnesses are equal.
    
    num_sources = size(xy0,1);              % K sources total
    xy = 2*sigma*(rand(num_sources,2)-.5);    % initial source positions
   
    % initialize expectation maximization parameters
    converged = 0;                               % convergence flag
    t = 1;                                       % iteration count
    max_iters = 5e2;                             % maximum gradient iterations
    
    xy_track = zeros(num_sources,2,max_iters+1);
    xy_track(:,:,t) = xy;
    while ~ converged
        % compute sample assignment weights
        weights = DirectImagingMeasurementProb(xy_samples,xy,sigma);
        weights = weights ./ sum(weights,2);

        % update the source position estimates
        xy = sum(permute(weights,[2,3,1]).* permute(xy_samples,[3,2,1]), 3)./sum(permute(weights,[2,3,1]),3);

        % update counter
        t = t+1;

        % add current estimate to track
        xy_track(:,:,t) = xy;
        
        % update convergence requirement
        converged = (t==max_iters) || all(vecnorm(xy_track(:,:,t)-xy_track(:,:,t-1),2,2) < sigma*1e-4);

    end
end

function p = DirectImagingMeasurementProb(xy,xy0,sigma)
    % returns the probability of detecting a photon at the location xy
    % given it was emitted by a source at position xy0.
    
    num_sources = size(xy0,1);
    num_samples = size(xy,1);

    p = zeros(num_samples,num_sources);
    for k = 1:num_sources
        p(:,k) = mvnpdf(xy,xy0(k,:),eye(2)*sigma^2);
    end
end