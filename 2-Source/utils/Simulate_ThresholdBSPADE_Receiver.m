function [xsk_est, outdata] = Simulate_ThresholdBSPADE_Receiver(xsk_0,N,sigma,dm,target_std_s)
% Description:
% Simulate a 2-stage receiver where one allocates as many resources as
% needed to the calibration stage in order to guarantee good localization.

% extract 2-source parameters
x0 = xsk_0(1);
s0 = xsk_0(2);
k0 = xsk_0(3);


%% STAGE 1: CALIBRATION
% Simulate a direct imaging measurement (equal brightness)

% setup switching criteria
M1 = 0;
x = [];
mom2_est = 0;
mom2_est_sig = inf;
dm_step = dm;

% while the second moment estimate of the direct imaging 
% distribution minus 1 is less than the standard deviation of the
% second moment estimate, continue integrating photons.
while (mom2_est - 1) < mom2_est_sig
    % collect photons
    x_next = DirectImagingMeasurement(x0,s0,0,dm_step,sigma);
    x = [x; x_next];
    M1 = numel(x);

    % Estimate midpoint and variance (centered second mom)
    x0_est = mean(x);
    mom2_est = mean(((x-x0_est)/sigma).^2);
    mom2_est_sig = sqrt(2*(M1-1)/(M1^2));
    
    % update the step size to meet the variance requirement
    tau = 2/(mom2_est-1)^2;
    discriminant= sqrt(tau^2 - 4*tau);
    M1_req = tau/2 + [discriminant,-discriminant]/2;
    M1_req = ceil(max(M1_req));
    dm_step = max(M1_req-M1,dm);
end

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

% make domain for centroid parameter
e = sigma/(sqrt(M1))*linspace(-3,3,71)';           de = e(2)-e(1);
e = permute(e,[2,1]);

% centroid prior distribution
p_e = normpdf(e,0,sigma/sqrt(M1));

% Simulate a binary spade measurement (equal brightness)
M2 = 0;
q = 0;
initial = 1;
while (initial) || (std_s > target_std_s)

    % simulate a BSPADE measurement
    q_next = BSPADEMeasurement(x0-x0_est,s0,0,dm,sigma);
    q = q + q_next;
    M2 = M2 + dm;

    % make domain for separation parameter
    min_s = max(mean_s-4*std_s,1e-9); max_s = mean_s+4*std_s;
    s = linspace(min_s,max_s,71)';  ds = s(2)-s(1);
    s = permute(s,[3,2,1]);

    % separation prior distribution
    p_s = gampdf((s/2/sigma).^2,alpha,1/beta).*(s/2/(sigma^2));

    % compute BSPADE likelihood
    BSPADE_likelihood = BSPADELikelihood(q,M2,e,s,0,sigma,1);

    % get the conditional posterior distribution on the separation P(s|epsilon,q)
    P_sIe = BSPADE_likelihood.*p_s;
    P_sIe = P_sIe./(sum(P_sIe,3)*ds + 1e-20);
        
    % marginal distribution on separation
    P_s = sum(P_sIe.*p_e*de,2);

    % calculate updated standard deviation estimate of separation under
    % marginal posterior
    mean_s = sum(P_s.*s,3)*ds;
    mean_s2 = sum(s.^2.*P_s,3)*ds;
    var_s = mean_s2-mean_s^2;
    std_s = sqrt(var_s);
    
    % no longer initial
    initial = 0;
end

%% STAGE 2: ESTIMATION
% make domain brightness parameter
k = linspace(-.5,.5,71)';   dk = k(2)-k(1);
k = permute(k,[4,3,2,1]);

% simulate a direct imaging measurmement (biased brightness)
xx = DirectImagingMeasurement(x0-x0_est,s0,k0,N,sigma);

% compute direct imaging likelihood
DI_likelihood = DirectImagingLikelihood_approx(xx,e,s,k,sigma,1);
DI_likelihood(isinf(DI_likelihood)) = max(DI_likelihood(~isinf(DI_likelihood)),[],'all'); % remove infs

% get the conditional posterior distribution on the brightness
P_kIse = DI_likelihood./(sum(DI_likelihood,4)*dk + 1e-20);

% get the marginal distribution on the brightness
P_k = sum(P_kIse.*P_sIe.*p_e*ds*de,[2,3]);

% Get MMSE estimator for separation
s_mmse = sum(P_s.*s*ds,3);

% Get MMSE estimator for brightness
k_mmse = sum(P_k.*k*dk,4);

% get conditional Max likelihood estimator for brightness
k_mle = mean(xx)/(2*s_mmse);
k_mle = min(max(-.5,k_mle),.5);

% collect the estimates
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
end