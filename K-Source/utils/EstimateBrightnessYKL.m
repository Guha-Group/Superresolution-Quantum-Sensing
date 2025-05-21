function b_est = EstimateBrightnessYKL(n,YKL,Psi_est)
% n         :   [K+1,1] vector of photon counts per YKL mode
% YKL       :   [K,K] matrix of projectors onto YKL measurements.
% Psi_est   :   [K,K] matrix of pure state vectors (each column is a unique
%               vector) corresponding to the field supplied by each point
%               source.


% number of states
num_states = size(Psi_est,1);

% conditional probability matrix providing the probability of detection in
% each YKL mode given emission by a particular source.
Q = abs(YKL'*Psi_est).^2;

% partition measured photons
ns = n(1:end-1); % support photons
n0 = n(end);     % bucket photons
ps = ns/sum(ns); % target distribution

% setup manifold optimization
D_KL = @(p,q) p'*log(p./q); % Kullback-Leibler Divergence
problem.M = multinomialfactory(num_states,1);   % manifold is the standard simplex
problem.cost = @(x) D_KL(ps,Q*x);               % cost is the KL divergence
%problem.egrad = @(x) Q.'*(ns./(Q*x));          % euclidean gradient of the KL cost

% solve for the brightness estimate
options.verbosity = 0;
warning('off', 'manopt:getGradient:approx');
warning('off', 'manopt:getHessian:approx');
[b_est,~,~,~] = trustregions(problem,[],options);
end