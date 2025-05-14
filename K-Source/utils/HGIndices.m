function [n,m] = HGIndices(n_max)
% returns the collection of 2D HG mode indices up to the level where
% n+m <= n_max.
n = zeros((n_max+1)*(n_max+2)/2,1);   % x-axis indices
m = n;                                % y-axis indices
k=1;
for i = 0:n_max
    for j = 0:i
        n(k) = j; m(k)=i-j;
        k=k+1;
    end
end
end