function [xsk_est,PDF] = SimulateDirectImagingReceiver(xsk_0,sigma,M,N)
% Description: Simulates a sequential estimation protocol for
% determining the centroid, separation, and brightness of two
% incoherently radiating sub-diffraction point sources. The receiver runs 
% as follows.
% 
%% Stage 1: Calibration
% M photons are detected from a scene with equal source brightnesses. The
% receiver constructs a prior distribution for the centroid and separation
% parameter under the assumption that the scene is sub-diffraction. 
%
%% Stage 2: Estimation
% N photons are detected from a scene with unbalanced source brightnesses.
% Given the prior distributions for the centroid and separation generated 
% in the first stage, we compute the marginal posterior for the 
% brightness. Ultimately, the receiver returns the Maximum likelihood 
% estimate of the centroid, the minimum mean squared error estimator of the
% separation, and the mean of the marginal posterior for the brightness parameter.
%
%%%%%%%%%%% INPUTS %%%%%%%%%%%%%
%   xsk_0   :   [1,3] vector of ground truth parameters [x0,s0,k0]
%   M       :   Number of photons allocated to Stage 1. If the
%               adaptive flag is false, all M photons are allocated to
%               the BSPADE measurement and we integrate as many photons
%               as needed in the direct imaging phase of the calibration 
%               to guarantee good performance with the Binary SPADE 
%               stage. If the adaptive flag is true, then the M photons
%               are adaptively allocated between Direct Imaging and BSPADE
%   N       :   Number of photons allocated to Stage 2.
%
%%%%%%%%%%% OUTPUTS %%%%%%%%%%%%
%   xsk_est :   [1,3] vector of parameter estimates [x0_mle,s0_mmmse,k0_mmse]
%   PDF     :   A structure with all priors and marginalized posterior
%               distributions of the estimators and their domains.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author(s):    Nico Deshler, ndeshler@arizona.edu
% institution:  University of Arizona, Wyant College of Optical Sciences
% date:         April 1, 2025

% extract 2-source parameters
x0 = xsk_0(1);
s0 = xsk_0(2);
k0 = xsk_0(3);

%% STAGE 1: CALIBRATION
% Simulate a direct imaging measurement (equal brightness)
x = SimulateDirectImaging(x0,s0,0,M,sigma);

% Estimate midpoint and variance (centered second mom)
x0_est = mean(x);
mom2_est = mean(((x-x0_est)/sigma).^2);
mu = max(mom2_est-1,1e-9);

% Estimate separation prior hyperparameters
alpha =  mu^2* (M^2/(M-1)) /2;
beta = 2*mu  * (M^2/(M-1));
mean_s = 2*sigma*gamma(1/2+alpha)/(sqrt(beta)*gamma(alpha));
mean_s2 = 4*sigma^2*alpha/beta;
var_s = mean_s2 - mean_s^2;

%% STAGE 2: ESTIMATION
% make domains for all three parameters.
min_s = max(mean_s-4*sqrt(var_s),1e-9); max_s = mean_s+4*sqrt(var_s);
e = sigma/(sqrt(M))*linspace(-3,3,71)';           de = e(2)-e(1);
s = linspace(min_s,max_s,71)';                    ds = s(2)-s(1);
k = linspace(-.5,.5,71)';                         dk = k(2)-k(1);

% make each variable in a new array dimension
e = permute(e,[2,1]);
s = permute(s,[3,2,1]);
k = permute(k,[4,3,2,1]);

% simulate a direct imaging measurmement (biased brightness)
xx = SimulateDirectImaging(x0,s0,k0,N,sigma);

% compute direct imaging likelihood
DI_likelihood = DirectImagingLikelihood_approx(xx,e,s,k,sigma,1);
DI_likelihood(isinf(DI_likelihood)) = max(DI_likelihood(~isinf(DI_likelihood)),[],'all'); % remove infs

% prior distributions
p_e = normpdf(e,0,sigma/sqrt(M));
p_s = gampdf((s/2/sigma).^2,alpha,1/beta).*(s/2/(sigma^2));

% get the conditional posterior distribution on the brightness
P_kIse = DI_likelihood./(sum(DI_likelihood,4)*dk);

% get the marginal distribution on the brightness
P_k = sum(P_kIse.*p_s.*p_e*ds*de,[2,3]);

% Get MMSE estimator for separation
s_mmse = mean_s;

% Get MMSE estimator for brightness
k_mmse = sum(P_k.*k*dk,4);

% collect outputs estimates
xsk_est = [x0_est,s_mmse,k_mmse];

% and distributions
PDF.xsk_0 = xsk_0;         % ground truth params
PDF.xsk_est = xsk_est;     % estimated params
PDF.k0 = k0;         % ground truth k
PDF.x = e+x0_est;    % domain of x0 estimator
PDF.s = s;           % domain of s estimator
PDF.k = k;           % domian of k estimator
PDF.P_x0 = p_e;      % posterior distribution on x0 estimator
PDF.P_s = p_s;       % prior on s estimator
PDF.P_k = P_k;       % marginal posterior distribution on k estimator
PDF.photons = [M,N]; % photons allocated to each stage

end