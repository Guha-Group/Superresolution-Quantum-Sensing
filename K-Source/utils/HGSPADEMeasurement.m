function mode_counts = HGSPADEMeasurement(xyb,nm,n,sigma,rot_angles)
% Returns a random instance of an HG SPADE measurement consisting of n 
% photons under the scene given by xyb.

% get the mode detection probabilities
[p_nm,~] = HermiteGaussMeasurementProb(xyb,nm,sigma,rot_angles);

% sample from the mode detection probabilities
mode_counts = zeros(size(p_nm));
T = size(p_nm,4);
for t = 1:T
    mode_counts(1,1,:,t) = mnrnd(floor(n/T),squeeze(p_nm(1,1,:,t)));
end
% add remainder photons into the last rotation
mode_counts(1,1,:,t) = mode_counts(1,1,:,t) + permute(mnrnd(rem(n,T),squeeze(p_nm(1,1,:,t))),[1,3,2]);
end

