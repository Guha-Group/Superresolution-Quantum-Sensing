function err = LocalizationError(xy_src,xy_est)
    % Computes the average localization error per source for a given
    % estimate of the constellation
    
    % calculate the total error distance between all possible pairings of
    % the ground truth sources and the estimated sources
    num_sources = size(xy_src,1);
    P = perms(1:num_sources);   % matrix of all permutations of the source indices
    
    sum_dist = zeros(size(P,1),1);
    
    for j = 1:size(P,1)
         index_permutation = P(j,:);
         xy_delta = xy_src - xy_est(index_permutation,:);
         deltas = vecnorm(xy_delta,2,2); % a vector of Euclidean distances between each source and its estimate
         sum_dist(j) = sum(deltas);
    end
       
    
    % select the minimum cumulative distance as the true error.
    % this assumes that each estimated source location is matched to the
    % nearest ground truth source
    err = min(sum_dist);
       
    % get the average error per source by dividing the number of sources.
    err = err / num_sources;
    
end

