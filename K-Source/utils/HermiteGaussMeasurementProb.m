function [p_nm,grad_p_nm] = HermiteGaussMeasurementProb(xyb,nm,sigma,rot_angles)
% Description:
% Evaluates the probability of photon detection in the Hermite-Gauss mode
% basis under a constellation of point sources assuming an imaging system
% with a gaussian PSF. We allow for rotated HG measurements 
% so as to eliminate degenerate solutions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    INPUTS  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xyb           :   [K,3] array of positions and relative brightnesses of the
%                   points sources in 2D. Each row is (x_i,y_i,b_i)
% nm            :   [M,2] array of HG mode indices
% sigma         :   Characteristic width of the gaussian modes
% rot_angles    :   [T,1] array of scene rotation angles to cover (in radians)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   OUTPUTS  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p_nm          :   [1,1,M,T]   probability distribution over the modes for both
%                   measurements
% grad_p_nm     :   [K,3,M,T] gradient of the probability distribution with
%                   respect to the scene parameters. In particular:
%                       grad_p_nm(i,1,:,:) = d/d_x_i p_nm(1,1,:,:)
%                       grad_p_nm(i,2,:,:) = d/d_y_i p_nm(1,1,:,:)
%                       grad_p_nm(i,3,:,:) = d/d_b_i p_nm(1,1,:,:)

% extract parameters
b = xyb(:,3);
xy = xyb(:,1:2);
    
% extract and shape format mode indices
n = permute(nm(:,1),[3,2,1]);
m = permute(nm(:,2),[3,2,1]); 

% define HG mode prob function and its derivative
p = @(u,k) exp(-u).*(u.^k)./factorial(k);  % 1D poisson distr handle
dp = @(u,k) p(u,k).*(k./u - 1);            % derivative of 1D poiss distr.

% define rotation matrix
R = @(theta) [cos(theta),sin(theta); -sin(theta),cos(theta)];

% make containers for outputs
num_sources = size(xyb,1);          % K
num_modes = numel(n);               % M
num_rotations = numel(rot_angles);  % T
p_nm = zeros(1,1,num_modes,num_rotations);
grad_p_nm = zeros(num_sources,3,num_modes,num_rotations);

for t = 1:num_rotations

    % get rotated scene coordinates
    xy_t = xy*R(rot_angles(t));

    % get HG poisson coeffs associated with each source position
    ux = (xy_t(:,1)/(2*sigma)).^2;
    uy = (xy_t(:,2)/(2*sigma)).^2;

    % compute the HG mode probability
    p_nm(1,1,:,t) = sum(b.*p(ux,n).*p(uy,m),1);
    
    % compute the gradient of the mode probabilities with respect to each param (x_i, y_i , b_i) 
    grad_p_nm(:,1,:,t) = b.*dp(ux,n).*p(uy,m).*(xy_t(:,1)/(2*sigma^2));  % dp_nm / dx_i
    grad_p_nm(:,2,:,t) = b.*p(ux,n).*dp(uy,m).*(xy_t(:,2)/(2*sigma^2));  % dp_nm / dy_i
    grad_p_nm(:,3,:,t) = p(ux,n).*p(uy,m);      % dp_nm / db_i

end
end
