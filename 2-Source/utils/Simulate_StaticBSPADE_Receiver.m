function [xsk_est, outdata] = Simulate_StaticBSPADE_Receiver(xsk_0,N,M,sigma,splitting_ratio)
% Description: Simulates measurements for a Bayesian receiver that
% estimates the centroid, separation, and brightness of two radiating 
% sub-diffraction point sources. The receiver is static in that the user
% defines the 
% The algorithm runs as follows:
%
%%  Stage 1: Calibration
%   This stage is broken into two sequential phases (1) Direct Imaging,
%   (2) Binary SPADE. The first phase is dedicated to estimating the  
%   geometric midpoint x0 between the two sources and performing a
%   preliminary estimate of their separation. The second phase is 
%   dedicated to refining the estimating of the separation with a binary
%   SPADE measurement. Throughout the entire calibration stage, the 
%   sources are configured to be equally bright.
%
%%  Stage 2: Estimation
%   In this stage, the receiver reverts to a direct imaging
%   measurement. The sources assume an unequal brightness. The receiver
%   looks to estimate the brightness bias k0 conditioned on the
%   statistics of the estimators for the midpoint x0 and separation s0.
%
%%%%%%%%%%% INPUTS %%%%%%%%%%%%%
%   xsk_0   :   [1,3] vector of ground truth parameters [x0,s0,k0]
%   N       :   Number of photons allocated to Sensing Stage.
%   M       :   Number of photons allocated to Calibration Stage. If the
%               adaptive flag is false, all M photons are allocated to
%               the BSPADE measurement and we integrate as many photons
%               as needed in the direct imaging phase of the calibration 
%               to guarantee good performance with the Binary SPADE 
%               stage. If the adaptive flag is true, then the M photons
%               are adaptively allocated between Direct Imaging and BSPADE
%   splitting_ratio   :   fraction of calibration stage photons allocated to direct 
%               imaging v. BSPADE where 0<alpha<=1 determines the fraction
%               of photons allocated direct imaging. If alpha=1, this is
%               equivalent to a fully direct imaging receiver.
%%%%%%%%%%% OUTPUTS %%%%%%%%%%%%
%   xsk_est :   [1,3] vector of parameter estimates [x0_mle,s0_mmmse,k0_mmse]
%   outdata     :   A structure with all priors and marginalized posterior
%               distributions of the estimators and their domains.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author(s):    Nico Deshler, ndeshler@arizona.edu
% institution:  University of Arizona, Wyant College of Optical Sciences
% date:         March 27, 2025

% extract 2-source parameters
x0 = xsk_0(1);
s0 = xsk_0(2);
k0 = xsk_0(3);

% get number of photons allocated to direct imaging and BSPADE in
% calibration stage
M1 = ceil(splitting_ratio*M);
M2 = M-M1;

%% STAGE 1: CALIBRATION
% Simulate a direct imaging measurement (equal brightness)
x = DirectImagingMeasurement(x0,s0,0,M1,sigma);

% Estimate midpoint and variance (centered second mom)
[x0_est,mom2_est] = DirectImagingParameterPreEstimates(x,sigma);

% Estimate separation prior hyperparameters
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

% make domains for all three parameters.
samp_dim = 81; % samples per parameter
min_s = max(mean_s-4*std_s,1e-9); 
max_s = mean_s+4*std_s;
e = sigma/(sqrt(M1))*linspace(-3,3,samp_dim)';          de = e(2)-e(1);
s = linspace(min_s,max_s,samp_dim)';                    ds = s(2)-s(1);
k = linspace(-.5,.5,samp_dim)';                         dk = k(2)-k(1);

% make each variable in a new array dimension
e = permute(e,[2,1]);
s = permute(s,[3,2,1]);
k = permute(k,[4,3,2,1]);

% if photons are available for BSPADE
if M2>0
    % Simulate a binary spade measurement (equal brightness)    
    q = BSPADEMeasurement(x0-x0_est,s0,0,M2,sigma);

    % compute BSPADE likelihood
    BSPADE_likelihood = BSPADELikelihood(q,M2,e,s,0,sigma,1);
end

%% STAGE 2: ESTIMATION

% prior distributions
p_e = normpdf(e,0,sigma/sqrt(M1));
p_s = gampdf((s/2/sigma).^2,alpha,1/beta).*(s/2/(sigma^2));

% get the conditional posterior distribution on the separation
% P(s|epsilon,q) if BSPADE measurements were made
if M2>0
    P_sIe = BSPADE_likelihood.*p_s;
    P_sIe = P_sIe./(sum(P_sIe,3)*ds + 1e-20);
else
    P_sIe = p_s;
end

% marginal distribution on separation
P_s = sum(P_sIe.*p_e,2)*de;
P_s = P_s./(sum(P_s)*ds); % normalize

% simulate a direct imaging measurement (biased brightness)
xx = DirectImagingMeasurement(x0-x0_est,s0,k0,N,sigma);

%{ 
%ORIGINAL
% compute direct imaging likelihood
DI_likelihood = DirectImagingLikelihood_approx(xx,e,s,k,sigma,1);
DI_likelihood(isinf(DI_likelihood)) = max(DI_likelihood(~isinf(DI_likelihood)),[],'all'); % remove infs

% get the conditional posterior distribution on the brightness
P_kIse = DI_likelihood./(sum(DI_likelihood,4)*dk + 1e-20);

% get the marginal distribution on the brightness
P_k = sum(P_kIse.*P_sIe.*p_e,[2,3])*ds*de;
P_k = P_k./(sum(P_k)*dk); % normalize
%}

%DI_likelihood = DirectImagingLikelihood_approx(xx,e,s,k,sigma,0);
%DI_likelihood(isinf(DI_likelihood)) = max(DI_likelihood(~isinf(DI_likelihood)),[],'all'); % remove infs

%LL = sum(DI_likelihood.*P_sIe.*p_e,[2,3])*ds*de;
%P_k = LL;

%%
%{
S1=0;
S2=0;
for i = 1:(numel(xx)-1)
    S2 = S2 + sum(xx(i) .* xx((i+1):end));
    S1 = S1 + sum(xx(i) +  xx((i+1):end));
end

DI_likelihood_approx = 1 + (2.*k.*s/ (sigma)^2).*(sum(xx,1)-M*e) + (2.*k.*s/ (sigma)^2).^2 .* (S2 + S1 * e + (M*(M-1)/2) * e.^2);
P_kIse = DI_likelihood_approx./(sum(DI_likelihood_approx,4)*dk + 1e-20);


% approximate conditional posterior on the brightness
%xi = 2*s.*(sum(xx,1) - numel(xx)*e)/sigma^2;
%P_kIse = exp(xi .* k)./((2./xi) .* sinh(xi/2));

% conditional mmse estimator
%k_mmse_cond = -1/xi + (1/2)*coth(xi/2);
%k_mmse_avg = sum(k_mmse_cond.*P_sIe.*p_e,[2,3])*ds*de;
%}


% Get MMSE estimator for separation
s_mmse = sum(P_s.*s,3)*ds;

% Get MMSE estimator for brightness
%k_mmse = sum(P_k.*k,4)*dk;

% get conditional Max likelihood estimator for brightness
k_mle = mean(xx)/(2*s_mmse);
k_mle = min(max(-.5,k_mle),.5);
P_k = ones(size(k));

% get map for integrated likelihood under uniform
%LL = sum(DI_likelihood.*P_sIe.*p_e,[2,3])*ds*de;
%[~,map_id] = max(LL);
%k_map = k(map_id);

% collect estimates
xsk_est = [x0_est,s_mmse,k_mle];

% and distributions
outdata.xsk_0 = xsk_0;         % ground truth params
outdata.xsk_est = xsk_est;     % estimated params
outdata.x = e+x0_est;          % domain of x0 estimator
outdata.s = s;                 % domain of s estimator
outdata.k = k;                 % domian of k estimator
outdata.P_x0 = p_e;            % posterior distribution on x0 estimator
outdata.P_s = P_s;             % marginal distribution on s estimator
outdata.P_k = P_k;             % marginal distribution on k estimator
outdata.photons = [M1,M2,N];   % photons allocated to each stage
if M2>0
    outdata.p_s = p_s;         % prior on s estimator
end

end