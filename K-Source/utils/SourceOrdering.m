function [xy_out,d_out] = SourceOrdering(xy,xy0)
% orders the permutation of the estimated sources such that they minimize
% the cumulative pair-wise Euclidean distance with from the ground truth.

% compute cumulative euclidean distance for all orderings
ordering = perms(1:size(xy,1));
d = zeros(size(ordering,1),1);
for j = 1:size(ordering,1)
    d(j) = sum(vecnorm(xy(ordering(j,:),:) - xy0,2,2),1);
end

% get the minimum distance ordering
[d_out,j_out] = min(d); 
xy_out = xy(ordering(j_out,:),:);
end