function xy = HGLocalizeSources(xy0,nm,mode_counts,sigma,rot_angles)
    % Description:
    % Estimates the position of multiple point sources in a 2D field of view
    % from the measurement outcomes mode_counts and measurement mode
    % p_model. The sources are assumed to be equally bright. We run a
    % maximum likelihood estimator via gradient ascent.
    % 
    %%%%%%%%%% INPUTS %%%%%%%%%%%%
    % xy0           :   [K,2] matrix of ground-truth cartesian coordinates
    %                   for each source (use for degeneracy breaking)
    % mode_counts   :   [1,1,M,D] a vector containing the number of photons
    %                   measured in each of the M modes for D unique povms. 
    % p_model       :   @fn(xy) function handle to the measurement
    %                   probabilities. Returns an array of size [1,1,M,D]
    % grad_p_model  :   @fn(xy) function handle to the gradient of the mode
    %                   measurement probability w.r.t each source
    %                   coordinate. Returns an array of size [2,K,M,D]
    %%%%%%%%%% OUTPUTS %%%%%%%%%%%%
    % xy            :   [K,2] matrix of estimated cartesian coordinates for
    %                   all of the sources
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author(s): Nico Deshler
    % Date: April 8, 2025

    % initialize source position estimates
    num_sources = size(xy0,1);              % K sources total
    xy = sigma*(rand(num_sources,2)-.5);    % initial source positions
    b = ones(num_sources,1)/num_sources;    % initial source brightnesses
    xy = xy - sum(b.*xy,1);                 % align to center of intensity

    % initialize gradient ascent parameters
    converged = 0;                               % convergence flag
    t = 1;                                       % iteration count
    %learning_rate = 1/sum(mode_counts(:));       % ascent rate
    learning_rate = 1e-4;       % ascent rate (step size)
    max_iters = 5e4;                             % maximum gradient iterations
    
    % tracking objects
    loglike_track = zeros(1,max_iters);         % track the log likelihood
    len_gradvec_track = zeros(1,max_iters);     % track the magnitude of the log likelihood gradient vector
    xy_track = zeros(num_sources,2,max_iters+1);% track of the estimates

    % run gradient ascent with ADAM optimization
    beta1 = .9; 
    beta2 = .999;
    eps = 1e-10;
    mt = 0; 
    vt = 0; 
   
    % normalize the mode counts 
    pm = mode_counts./sum(mode_counts,3);

    while ~converged
        
        % compute the mode probabilities and gradients
        [p,grad_p] = HermiteGaussMeasurementProb([xy,b],nm,sigma,rot_angles);
        grad_xy_p = grad_p(:,1:2,:,:);

        % compute the log likelihood
        loglike = sum(pm.*log(p),[3,4]); % + log(multinomial coeffient);

        % compute the gradient of the log likelihood
        grad_xy_loglike = sum(pm.*grad_xy_p./p ,[3,4]);

        % adam optimizer parameter updats
        mt = beta1*mt + (1-beta1) *grad_xy_loglike;
        vt = (beta2)*vt + (1-beta2)*(grad_xy_loglike).^2;
        mt_hat = mt/(1-beta1^t) ; vt_hat = vt/(1-beta2^t);

        % update estimates
        xy = xy + mt_hat./(sqrt(vt_hat)+eps)*learning_rate;
        %xy = xy + grad_xy_loglike*learning_rate; % vanilla gradient descent

        % update convergence requirement
        len_gradvec = sqrt(grad_xy_loglike(:)'*grad_xy_loglike(:));
        converged = (len_gradvec < 1e-6) || (t==max_iters);


        %% update tracks
        loglike_track(t) = loglike;
        len_gradvec_track(t) = len_gradvec;
        xy_track(:,:,t) = xy;
        
        % update iteration counter
        if ~converged
            t = t + 1;
        end
    end

    % break degeneracies in final estimate
    xy = BreakDegeneracies(xy,xy0);
  
        
% -------------------%
% debugging plot
figure
tiledlayout(1,3,'Padding','compact','TileSpacing','compact')
nexttile(1)
plot(1:t,loglike_track(1:t),'k','LineWidth',1.5)
axis square; box on; grid on;
xlabel('iter'); ylabel('Log Likelihood');
set(gca,'yscale','log'); set(gca,'xscale','log'); 
title('Log-Likelihood Trace','Interpreter','latex');

nexttile(2)
plot(1:t,len_gradvec_track(1:t),'k','LineWidth',1.5)
xlabel('iter'); ylabel('Gradient Vector Length'); axis square; box on; grid on;
set(gca,'yscale','log'); %set(gca,'xscale','log'); 
title('Gradient Vector Length Trace','Interpreter','latex');

nexttile(3)
id = 1:100:t;
constellation_colors = turbo(numel(id));
hold on
for k = 1:numel(id)
    scatter(xy_track(:,1,id(k)),xy_track(:,2,id(k)),20,constellation_colors(k,:),...
        'filled','MarkerEdgeAlpha',1,'MarkerFaceAlpha',exp(-5* (1- (k/numel(id)))));
end
scatter(xy(:,1),xy(:,2),20,'r','d','filled')
scatter(xy0(:,1),xy0(:,2),20,'k','filled')
hold off

axis square; box on; grid on;
xticks(-1:.5:1); yticks(-1:.5:1);
xlim([-1,1]); ylim([-1,1]);
xlabel('$x/\sigma$','interpreter','latex')
ylabel('$y/\sigma$','interpreter','latex')
title('Source Location Estimates','interpreter','latex')
%--------------------%
end


function xy = BreakDegeneracies(xy,xy0)

    % compute all inversion symmetries of the source coordinates
    xy_set = GenerateDegeneracies(xy);
    
    %{
    % pick the one that has center of intensity closest to optical axis
    centers = mean(xy_set,1);
    d0 = vecnorm(centers,2,2);
    xy_set = xy_set(:,:,d0==min(d0));
    %}

    % then pick the the constellation that is closest to the ground truth
    label = perms(1:size(xy,1));
    d1 = zeros(1,size(xy_set,3));
    for i = 1:size(xy_set,3)
        temp = zeros(1,size(label,1));
        for j = 1:size(label,1)
            % compute the euclidean distance between the i'th constellation
            % with labelling permutation label(j,:) and the ground truth.
            temp(j) = sum(vecnorm(xy_set(label(j,:),:,i) - xy0,2,2),1);
        end
        d1(i) = min(temp);
    end
    xy = xy_set(:,:,d1==min(d1));
end

function xy_set = GenerateDegeneracies(xy)
    % Description: Returns every degenerate solution of maximum likelihood
    % estimation with spatial mode sorting. In particular, each estimated
    % coordinate admits three other degeneracy points corresponding to
    % reflection about x-axis, reflection about y-axis and reflection about
    % x then y.
   
    num_sources = size(xy,1);

    % put all coordinates in the positive quadrant
    xy = abs(xy); 

    % get the sign combinations for remaining quadrants
    sign_mat = [1,1;-1,1;-1,-1;1,-1];
    
    % determine all possible constellation
    cindex = GenerateCombinations(1:4,num_sources);
    num_constellations = size(cindex,1);
    xy_set = zeros(num_sources,2,num_constellations);
    for i = 1:num_sources
        xy_set(i,:,:) = sign_mat(cindex(:,i),:).';
    end
    xy_set = xy.*xy_set;

end

function C = GenerateCombinations(v,k)
    % Description: generates all vectors of length k whose entries are
    % sampled from vector v.
    
    % Generate all possible combinations
    grids = cell(1, k);
    [grids{:}] = ndgrid(v);
    
    % Convert grid output into a matrix where each row is a combination
    C = cell2mat(cellfun(@(x) x(:), grids, 'UniformOutput', false));
end