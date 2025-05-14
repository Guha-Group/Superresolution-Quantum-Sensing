function xy_samples = DirectImagingMeasurement(xyb,sigma,n)
%DIRECTIMAGINGMEASUREMENT Simulates a direct imaging measurement for a
% constellation of point sources imaged under a gaussian PSF.
%%%%%% INPUTS %%%%%%
% xyb       : [K,3] array containing the (x,y,b) coordinate of each source
%              (positions in 2D cartesian space as well as relative brightness)
% sigma     : Rayleigh scale
% n         : number of photons to sampel
%%%%%% OUTPUTS %%%%%%
% xy_samples: [n,2] array of (x,y) coordinates for each simulated photon
%             arrival location

% sample a multinomial distribution to get the number of photons emitted
% from each source
nk = mnrnd(n,xyb(:,3)); 

% sample the locations from a shifted gaussian
xy_samples = zeros(n,2);
nf = 0;
for k = 1:size(xyb,1)
    ns = nf + 1;
    nf = nf + nk(k);
    xy_samples(ns:nf,:) = mvnrnd(xyb(k,1:2),eye(2)*sigma^2,nk(k));
    
end
end

