sigma = 1;                                  % Rayleigh Limit (for gaussian PSF)
num_sources = 7;                            % number of sources in scene

% Ground Truth Parameters
sources_xy = sigma*(rand([num_sources,2])-.5);          % source cartesian coordinates
sources_b = rand(num_sources,1); 
sources_b = sources_b/sum(sources_b);
%Psi = generateIncoherentSourceStateVectors(sources_xy,sources_b,sigma);

% Estimated Parameters
sources_xy_est = sources_xy + sigma*(rand([num_sources,2])-.5)/10;
sources_b_est = sources_b;
Psi_est = generateIncoherentSourceStateVectors(sources_xy_est,sources_b_est,sigma);

TestMatGradPerr(Psi_est,sources_b)

% calculate the overlap integrals between estimated and ground truth
% states. To do so, compute the pairwise difference vectors between all sources
deltas_xy = permute(sources_xy,[1,3,2])-permute(sources_xy_est,[3,1,2]);

% Then compute the Pseudo-Gram Matrix (assuming gaussian PSF) containing 
% the overlaps <psi_est_i|psi_j>
T = exp(-sum((deltas_xy/sigma).^2,3)/8);

% get the Helstrom measurement under the estimated parameters
YKL = getYKL(Psi_est,sources_b_est);

% express the YKL measurement in the (non-orthogonal) state space basis
G_est = Psi_est'*Psi_est;
M = YKL'*Psi_est;
A = conj(M*inv(G_est));

% evaluate the overlaps between the YKL measurement projectors 
% and the ground truth pure states
YKL_state_overlaps = A*(T.');

% Evaluate the probability in each YKL projector mode
mode_prob = abs(YKL_state_overlaps).^2 * sources_b;

% append the bucket mode
HM_prob = [mode_prob; 1-sum(mode_prob)];


%% FIGURES
figure
tiledlayout(1+num_sources,num_sources,'TileSpacing','compact','Padding','compact')

% plot the ground truth and estimated source positions
nexttile(1,[num_sources,num_sources])
hold on
scatter(sources_xy(:,1),sources_xy(:,2),30*num_sources*sources_b,'k','filled')
scatter(sources_xy_est(:,1),sources_xy_est(:,2),30*num_sources*sources_b_est,'r','d','filled')
hold on
legend({'Ground Truth','Estimate'})
axis square
box on
grid on
xlim([-1,1])
ylim([-1,1])
xticks(-1:.5:1)
yticks(-1:.5:1)

% plot the YKL modes
[X,Y] = meshgrid(sigma*linspace(-5,5,1001));
PSF_modes = 1/sqrt(2*pi*sigma^2)*...
            exp(- (1/(2*sigma)^2)*(...
            (X-permute(sources_xy_est(:,1),[3,2,1])).^2 +...
            (Y-permute(sources_xy_est(:,2),[3,2,1])).^2 ));
colormap(turbo)
for t = 1:num_sources
    nexttile(num_sources^2+t)
    mode = sum(permute(A(t,:),[3,1,2]).*PSF_modes,3);
    imagesc([min(X(:)),max(X(:))]/sigma,[min(Y(:)),max(Y(:))]/sigma,abs(mode).^2);
    axis square
    set(gca,'ydir','normal')
end



%% FUNCTIONS 
function Psi = generateIncoherentSourceStateVectors(sources_xy,priors,sigma)
    % returns the representation of the incoherent source states in the
    % eigenbasis of the density operator.
    % 
    % Psi       : matrix of states. Each column is a unique state.
    
    % compute the pairwise difference vectors between all sources
    deltas_xy = permute(sources_xy,[1,3,2])-permute(sources_xy,[3,1,2]);
    
    % make diagonal matrix of prior probabilities
    P = diag(priors(:));

    % Gram Matrix (assuming gaussian PSF)
    G = exp(-sum((deltas_xy/sigma).^2,3)/8);
  
    % Get spectral decomposition of Gram matrix
    [U,D] = eig(G);
    
    % Determine the S matrix
    A = sqrt(P)*U*sqrt(D);
    S = A'*A;
    
    % Get spectral decomposition of S matrix (the eigenvalues of S are the same
    % as the eigenvalues of rho)
    [V_dagger, rho_evals] = eig(S);
    V = V_dagger';
    
    % Compute the representation matrix R (and its inverse) of Psi in the 
    % eigenbasis of rho
    R = V * sqrt(D) * U';
    Psi = R; % the states as column vectors represented in the eigenbasis of the density operator
end


function YKL = getYKL(Psi,priors)
    
    % the number of states
    num_states = numel(priors);
    
    % setup manifold optimization
    problem.M = unitaryfactory(num_states,1);                    % manifold is the space of unitaries (orthnormal projectors)
    problem.cost = @(x) 1- abs(sum(conj(x).*Psi,1)).^2 * priors; % cost is the probability of error success
    problem.egrad = @(x) - (sum(conj(Psi).*x,1).*priors.').*Psi; % DO NOT FUCK WITH THIS -- TOOK YOU FOREVER TO DERIVE THE RIGHT EUCLIDEAN GRAD
    
    % solve the optimal measurement
    [YKL,~,~,~] = trustregions(problem);
    
end

function TestMatGradPerr(Psi,priors)

    % the number of states
    num_states = numel(priors);


    P_err = @(x) 1 - abs(sum(conj(x).*Psi,1)).^2 * priors;
    gradP_err = @(x) -(x .* abs(Psi).^2).*(priors.');

    % evaluate the probability of error at some unitary
    temp = randn(num_states) + 1i*randn(num_states);
    temp = temp+temp';
    [U,~] = eig(temp);

    % consider a small displacement along the matrix M
    t = 1e-8;
    M =randn(num_states) + 1i*randn(num_states);

    differential_p = (P_err(U+t*M) - P_err(U))/t;

    % evaluate whether this is equal to the projection with the gradient or
    % not
    inner_product = @(X,Y) real(trace(X*Y'));
    ip = 2*inner_product(gradP_err(U),M);
end

