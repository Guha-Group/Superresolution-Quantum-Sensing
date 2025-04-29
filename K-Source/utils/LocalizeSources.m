function xy = LocalizeSources(xy0,mode_counts,p_model,grad_p_model,sigma)
    % Description:
    % Estimates the position of multiple point sources in a 2D field of view
    % from the measurement outcomes mode_counts and measurement mode
    % p_model. The sources are assumed to be equally bright. We run a
    % maximum likelihood estimator via gradient ascent.
    % 
    %%%%%%%%%% INPUTS %%%%%%%%%%%%
    % xy0           :   [2,K] matrix of ground-truth cartesian coordinates
    %                   for each source (use for degeneracy breaking)
    % mode_counts   :   [1,1,M,D] a vector containing the number of photons
    %                   measured in each of the M modes for D unique povms. 
    % p_model       :   @fn(xy) function handle to the measurement
    %                   probabilities. Returns an array of size [1,1,M,D]
    % grad_p_model  :   @fn(xy) function handle to the gradient of the mode
    %                   measurement probability w.r.t each source
    %                   coordinate. Returns an array of size [2,K,M,D]
    %%%%%%%%%% OUTPUTS %%%%%%%%%%%%
    % xy            :   [2,K] matrix of estimated cartesian coordinates for
    %                   all of the sources
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author(s): Nico Deshler
    % Date: April 8, 2025

    % initialize source position estimates
    num_sources = size(xy0,1); % K sources total
    xy = sigma*rand(size(xy)); % initial source positions
    xy = xy - mean(xy,1); % align to center of intensity


    % initialize gradient ascent parameters
    converged = 0;              % convergence flag
    t = 1;                        % iteration count
    alpha;                      % ascent rate (step size)
    max_iters = 500;            % maximum gradient iterations
    
    % tracking objects
    %logLikelihood_track = zeros(1,max_iters);   % track the log likelihood
    %len_gradvec_track = zeros(1,max_iters);     % track the magnitude of the log likelihood gradient vector
    %xy_track = zeros(num_sources,2,max_iters+1);% track of the estimates
    %xy_track(:,:,t)=xy;

    % run gradient ascent
    while ~converged
        
        % compute the log likelihood
        log_like = sum(mode_counts.*log(p_model(xy))-p_model(xy));

        % compute the gradient of the log likelihood
        grad_log_like = sum( mode_counts./p_model(xy) .* grad_p_model(xy) - grad_p_model(xy));

        % update estimates
        xy = xy + grad_log_like*learning_rate;

        % update iteration counter
        t = t + 1;

        % update convergence requirement
        converged;
    end
end