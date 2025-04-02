function [xsk_est,PDF] = SimulateBayesianReceiver(xsk_0,sigma,M,N,dm,adaptive)
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
    %   M       :   Number of photons allocated to Stage 1. If the
    %               adaptive flag is false, all M photons are allocated to
    %               the BSPADE measurement and we integrate as many photons
    %               as needed in the direct imaging phase of the calibration 
    %               to guarantee good performance with the Binary SPADE 
    %               stage. If the adaptive flag is true, then the M photons
    %               are adaptively allocated between Direct Imaging and BSPADE
    %   N       :   Number of photons allocated to Stage 2.
    %   dm      :   Number of photons to integrate between adaptations in
    %               Calibration stage
    %   adaptive:   Flag to indicate whether adaptive estimation should be
    %               used or not.
    %%%%%%%%%%% OUTPUTS %%%%%%%%%%%%
    %   xsk_est :   [1,3] vector of parameter estimates [x0_mle,s0_mmmse,k0_mmse]
    %   PDF     :   A structure with all priors and marginalized posterior
    %               distributions of the estimators and their domains.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % author(s):    Nico Deshler, ndeshler@arizona.edu
    % institution:  University of Arizona, Wyant College of Optical Sciences
    % date:         March 27, 2025

    % extract 2-source parameters
    x0 = xsk_0(1);
    s0 = xsk_0(2);
    k0 = xsk_0(3);
    
    %% STAGE 1: CALIBRATION
    % Simulate a direct imaging measurement (equal brightness)

    if ~adaptive %% NON-ADAPTIVE (TYPE-I ESTIMATION)
        
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
            x_next = SimulateDirectImaging(x0,s0,0,dm_step,sigma);
            x = [x; x_next];
            M1 = numel(x);
        
            % Estimate midpoint and variance (centered second mom)
            x0_est = mean(x);
            mom2_est = mean(((x-x0_est)/sigma).^2);
            mom2_est_sig = sqrt(2*(M1-1)/(M1^2));
            
            %{
            % update the step size to meet the variance requirement
            tau = 2/(mom2_est-1)^2;
            discriminant= sqrt(tau^2 - 4*tau);
            M1_req = tau/2 + [discriminant,-discriminant]/2;
            M1_req = max(M1_req);
            dm_step = round(max(M1_req-M1,dm));
            %}
        end
        
        % Estimate separation prior hyperparameters
        alpha =  (mom2_est-1)^2* (M1^2/(M1-1)) /2;
        beta = 2*(mom2_est-1)  * (M1^2/(M1-1));
        mean_s = 2*sigma*gamma(1/2+alpha)/(sqrt(beta)*gamma(alpha));
        mean_s2 = 4*sigma^2*alpha/beta;
        var_s = mean_s2 - mean_s^2;
        
        % Simulate a binary spade measurement (equal brightness)
        M2 = M;
        q = SimulateBSPADE(x0-x0_est,s0,0,M2,sigma);

    else         %% ADAPTIVE (TYPE-II ESTIMATION)

        M1 = 0;
        x = [];
        expected_var_s = zeros(floor(M/dm)+1,1);
        switch_condition = 0;
        iter = 0;
        while ~switch_condition && M1<M
            % collect photons
            x_next = SimulateDirectImaging(x0,s0,0,dm,sigma);
            x = [x;x_next];
            M1 = numel(x);
        
            % Estimate midpoint and variance of direct imaging PDF
            x0_est = mean(x);
            mom2_est = mean(((x-x0_est)/sigma).^2);
            mom2_est_sig = sqrt(2*(M1-1)/(M1^2));
            mom2_est = max(mom2_est,1+mom2_est_sig); % preserves stability 
        
            % Estimate separation prior parameters
            alpha =  (mom2_est-1)^2* (M1^2/(M1-1)) /2;
            beta = 2*(mom2_est-1)  * (M1^2/(M1-1));
            mean_s = 2*sigma*gamma(1/2+alpha)/(sqrt(beta)*gamma(alpha));
            mean_s2 = 4*sigma^2*alpha/beta;
            var_s = mean_s2 - mean_s^2;
        
            % make domains for midpoint and separation params
            min_s = max(mean_s-4*sqrt(var_s),1e-9); max_s = mean_s+4*sqrt(var_s);
            e = sigma/(sqrt(M1))*linspace(-3,3,71)';
            s = linspace(min_s,max_s,121)';                  
        
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
        q = SimulateBSPADE(x0-x0_est,s0,0,M2,sigma);
    end

    %% STAGE 2: ESTIMATION

    % make domains for all three parameters.
    min_s = max(mean_s-4*sqrt(var_s),1e-9); max_s = mean_s+4*sqrt(var_s);
    e = sigma/(sqrt(M1))*linspace(-3,3,71)';           de = e(2)-e(1);
    s = linspace(min_s,max_s,71)';                    ds = s(2)-s(1);
    k = linspace(-.5,.5,71)';                         dk = k(2)-k(1);
    
    % make each variable in a new array dimension
    e = permute(e,[2,1]);
    s = permute(s,[3,2,1]);
    k = permute(k,[4,3,2,1]);
 
    % simulate a direct imaging measurmement (biased brightness)
    xx = SimulateDirectImaging(x0-x0_est,s0,k0,N,sigma);
                   
    % compute BSPADE likelihood
    BSPADE_likelihood = BSPADELikelihood(q,M2,e,s,0,sigma);
    
    % compute direct imaging likelihood
    DI_likelihood = DirectImagingLikelihood_approx(xx,e,s,k,sigma,1);
    DI_likelihood(isinf(DI_likelihood)) = max(DI_likelihood(~isinf(DI_likelihood)),[],'all'); % remove infs
    
    % prior distributions
    p_e = normpdf(e,0,sigma/sqrt(M1));
    p_s = gampdf((s/2/sigma).^2,alpha,1/beta).*(s/2/(sigma^2));
    
    % get the conditional posterior distribution on the separation P(s|epsilon,q)
    P_sIe = BSPADE_likelihood.*p_s;
    P_sIe = P_sIe./(sum(P_sIe,3)*ds);
    
    % marginal distribution on separation
    P_s = sum(P_sIe.*p_e*de,2);
    
    % get the conditional posterior distribution on the brightness
    P_kIse = DI_likelihood./(sum(DI_likelihood,4)*dk);
    
    % get the marginal distribution on the brightness
    P_k = sum(P_kIse.*P_sIe.*p_e*ds*de,[2,3]);
    
    % Get MMSE estimator for separation
    s_mmse = sum(P_s.*s*ds,3);
    
    % Get MMSE estimator for brightness
    k_mmse = sum(P_k.*k*dk,4);

    % collect outputs estimates
    xsk_est = [x0_est,s_mmse,k_mmse];

    % show figures
    %Figures

    % and distributions
    PDF.xsk_0 = xsk_0;         % ground truth params
    PDF.xsk_est = xsk_est;     % estimated params
    PDF.k0 = k0;         % ground truth k
    PDF.x = e+x0_est;    % domain of x0 estimator
    PDF.s = s;           % domain of s estimator
    PDF.k = k;           % domian of k estimator
    PDF.p_s = p_s;       % prior on s estimator
    PDF.P_x0 = p_e;      % posterior distribution on x0 estimator
    PDF.P_s = P_s;       % marginal posterior distribution on s estimator
    PDF.P_k = P_k;       % marginal posterior distribution on k estimator
    PDF.photons = [M1,M2,N]; % photons allocated to each stage
end