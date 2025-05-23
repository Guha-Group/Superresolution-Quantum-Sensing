function [xsk_est, outdata] =Simulate_AdaptiveBSPADE_Receiver(xsk_0,N,M,sigma,dm)
% Description: Simulates measurements for a Bayesian receiver that
% estimates the centroid, separation, and brightness of two radiating 
% sub-diffraction point sources. The algorithm runs as follows:
%
%%   Stage 1: Calibration
%   This stage is broken into two sequential phases (1) Direct Imaging,
%   (2) Binary SPADE. The first phase is dedicated to estimating the  
%   geometric midpoint x0 between the two sources and performing a
%   preliminary estimate of their separation. The second phase is 
%   dedicated to refining the estimating of the separation with a binary
%   SPADE measurement. Throughout the entire calibration stage, the 
%   sources are configured to be equally bright.
%
%%   Stage 2: Estimation
%   In this stage, the receiver reverts to a direct imaging
%   measurement. The sources assume an unequal brightness. The receiver
%   looks to estimate the brightness bias k0 conditioned on the
%   statistics of the estimators for the midpoint x0 and separation s0.
%
%%%%%%%%%%% INPUTS %%%%%%%%%%%%%
%   xsk_0   :   [1,3] vector of ground truth parameters [x0,s0,k0]
%   N       :   Number of photons allocated to Stage 2.
%   M       :   Number of photons allocated to Stage 1. If the
%               adaptive flag is false, all M photons are allocated to
%               the BSPADE measurement and we integrate as many photons
%               as needed in the direct imaging phase of the calibration 
%               to guarantee good performance with the Binary SPADE 
%               stage. If the adaptive flag is true, then the M photons
%               are adaptively allocated between Direct Imaging and BSPADE
%   dm      :   Number of photons to integrate between adaptations in
%               Calibration stage
%%%%%%%%%%% OUTPUTS %%%%%%%%%%%%
%   xsk_est :   [1,3] vector of parameter estimates [x0_mle,s0_mmmse,k0_mmse]
%   outdata :   A structure with all priors and marginalized posterior
%               distributions of the estimators and their domains.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author(s):    Nico Deshler, ndeshler@arizona.edu
% institution:  University of Arizona, Wyant College of Optical Sciences
% date:         March 27, 2025

% extract 2-source parameters
x0 = xsk_0(1);
s0 = xsk_0(2);
k0 = xsk_0(3);

% samples per parameter
samp_dim = 121; 

% Simulate a direct imaging measurement (equal brightness)
M1 = 0;
x = [];
expected_var_s = zeros(floor(M/dm)+1,1);
switch_condition = 0;
iter = 0;
while ~switch_condition && M1<M
    % collect photons
    x_next = DirectImagingMeasurement(x0,s0,0,dm,sigma);
    x = [x; x_next];
    M1 = numel(x);
    
    % Estimate midpoint and variance of direct imaging PDF  
    [x0_est,mom2_est] = DirectImagingParameterPreEstimates(x,sigma);
    
    % Estimate separation prior parameters
    alpha =  (mom2_est-1)^2* (M1^2/(M1-1)) /2;
    beta = 2*(mom2_est-1)  * (M1^2/(M1-1));
    % use asymptotic form for the ratio of gamma functions for large alpha
    if isinf(gamma(alpha))
        % https://www.johndcook.com/blog/2018/03/14/approximating-gamma-ratios/
        temp = alpha - 1/2;
        gamma_ratio_approx = (temp^2 + temp/2 + temp/8)^(1/4); 
        gamma_ratio = gamma_ratio_approx;
    else
        gamma_ratio = (gamma(alpha+1/2)/gamma(alpha));
    end
    mean_s = (2*sigma/sqrt(beta))*gamma_ratio;
    mean_s2 = 4*sigma^2*alpha/beta;
    var_s = mean_s2 - mean_s^2;
    std_s = sqrt(var_s);
    
    % make domains for midpoint and separation params
    min_s = max(mean_s-4*std_s,1e-9); 
    max_s = mean_s+4*std_s;
    e = sigma/(sqrt(M1))*linspace(-3,3,samp_dim)';    de=e(2)-e(1);
    s = linspace(min_s,max_s,samp_dim)';              ds=s(2)-s(1); 
    
    % make each variable in a new array dimension
    e = permute(e,[2,1]);
    s = permute(s,[3,2,1]);
    
    % prior distributions
    p_e = normpdf(e,0,sigma/sqrt(M1));
    p_s = gampdf((s/2/sigma).^2,alpha,1/beta).*(s/2/(sigma^2));
    
    % perform variance look ahead to see if the estimate of
    % separation parameter is decreasing. If not, switch to BSPADE.
    E_var_s = AdaptiveVarianceLookAhead(M,M1,p_s,p_e,s,e,sigma);
    iter = iter + 1;
    expected_var_s(iter) = E_var_s;
    if iter>1
        switch_condition = expected_var_s(iter)>expected_var_s(iter-1);
    end
end

% Simulate a binary spade measurement (equal brightness)
M2 = M - M1;
q = BSPADEMeasurement(x0-x0_est,s0,0,M2,sigma);

%% STAGE 2: ESTIMATION
% make domain for brightness parameter
k = linspace(-.5,.5,samp_dim)';                         dk = k(2)-k(1);

% put it into a new array dimension
k = permute(k,[4,3,2,1]);

% simulate a direct imaging measurmement (biased brightness)
xx = DirectImagingMeasurement(x0-x0_est,s0,k0,N,sigma);
               
% compute BSPADE likelihood
BSPADE_likelihood = BSPADELikelihood(q,M2,e,s,0,sigma,1);

% compute direct imaging likelihood
DI_likelihood = DirectImagingLikelihood_approx(xx,e,s,k,sigma,1);
DI_likelihood(isinf(DI_likelihood)) = max(DI_likelihood(~isinf(DI_likelihood)),[],'all'); % remove infs

% get the conditional posterior distribution on the separation P(s|epsilon,q)
P_sIe = BSPADE_likelihood.*p_s;
P_sIe = P_sIe./(sum(P_sIe,3)*ds + 1e-20);

% marginal distribution on separation
P_s = sum(P_sIe.*p_e,2)*de;

% get the conditional posterior distribution on the brightness
P_kIse = DI_likelihood./(sum(DI_likelihood,4)*dk + 1e-20);

% get the marginal distribution on the brightness
P_k = sum(P_kIse.*P_sIe.*p_e,[2,3])*ds*de;

% Get MMSE estimator for separation
s_mmse = sum(P_s.*s,3)*ds;

% Get MMSE estimator for brightness
k_mmse = sum(P_k.*k,4)*dk;

% get conditional Max likelihood estimator for brightness
k_mle = mean(xx)/(2*s_mmse);
k_mle = min(max(-.5,k_mle),.5);

% collect outputs estimates
xsk_est = [x0_est,s_mmse,k_mmse];

% and distributions
outdata.xsk_0 = xsk_0;         % ground truth params
outdata.xsk_est = xsk_est;     % estimated params
outdata.x = e+x0_est;    % domain of x0 estimator
outdata.s = s;           % domain of s estimator
outdata.k = k;           % domian of k estimator
outdata.p_s = p_s;       % prior on s estimator
outdata.P_x0 = p_e;      % posterior distribution on x0 estimator
outdata.P_s = P_s;       % marginal posterior distribution on s estimator
outdata.P_k = P_k;       % marginal posterior distribution on k estimator
outdata.photons = [M1,M2,N]; % photons allocated to each stage
outdata.dm = dm;
outdata.iter = iter;
outdata.expected_var_s = expected_var_s;
end