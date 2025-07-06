function mode_counts = HGSPADEMeasurement(xyb,nm,n,sigma,rot_angles)
% Returns a random instance of an HG SPADE measurement consisting of n 
% photons under the scene given by xyb.

% get the mode detection probabilities
[p_nm,~] = HermiteGaussMeasurementProb(xyb,nm,sigma,rot_angles);

% adjust probability vector to ensure the entries sum to 1 (stability for
% mnrnd)
p_nm(1,1,all(nm==0,2),:) = p_nm(1,1,all(nm==0,2),:) - (sum(p_nm,3)-1); 
p_nm(1,1,all(nm==0,2),:) = p_nm(1,1,all(nm==0,2),:) - (sum(p_nm,3)-1); 

% sample from the mode detection probabilities
mode_counts = zeros(size(p_nm));
T = numel(rot_angles);
for t = 1:T
    mode_counts(1,1,:,t) = mnrnd(floor(n/T),squeeze(p_nm(1,1,:,t)));
end
% add remainder photons into the last rotation
if rem(n,T)
    mode_counts(1,1,:,T) = mode_counts(1,1,:,T) + permute(mnrnd(rem(n,T),squeeze(p_nm(1,1,:,T))),[1,3,2]);
end

% set nan elements to zero
mode_counts(isnan(mode_counts))=0;
end

