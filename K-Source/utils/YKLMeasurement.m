function [mode_counts,YKL,Psi_est] = YKLMeasurement(xyb,xyb_est,sigma,n)
% Simulates an n-photon measurement for the YKL basis under an estimate of
% a constellation.

% get estimated pure states represented in the eigenbasis of the estimated
% density operator
Psi_est = generateIncoherentSourceStateVectors(xyb_est,sigma);

% calculate the overlap integrals between estimated and ground truth
% states (assuming gaussian PSF). This generates a Pseudo-Gram Matrix 
% containing the overlaps <psi_est_i|psi_j>
deltas_xy = permute(xyb(:,1:2),[1,3,2]) - permute(xyb_est(:,1:2),[3,1,2]);
T = exp(-sum((deltas_xy/sigma).^2,3)/8);

% get the YKL/Helstrom measurement under the estimated parameters
YKL = getYKL(Psi_est,xyb_est(:,3));

% express the YKL measurement in the (non-orthogonal) state space basis
G_est = Psi_est'*Psi_est;
M = YKL'*Psi_est;
A = conj(M*inv(G_est));

% evaluate the overlaps between the YKL measurement projectors 
% and the ground truth pure states
YKL_state_overlaps = A*(T.');

% Evaluate the probability in each YKL projector mode
mode_prob = abs(YKL_state_overlaps).^2 * xyb(:,3);

% append the bucket mode to collect residual probability arising from
% position estimate mismatch
mode_prob = [mode_prob; 1-sum(mode_prob)];

% sample from YKL distribution probability
mode_counts = mnrnd(n,mode_prob)';

%%%%%%%%%%%%%%%%%%%%%%%
%% Debugging Figures %%
%%%%%%%%%%%%%%%%%%%%%%%
figure
num_sources = size(xyb,1);
tiledlayout(1+num_sources,num_sources,'TileSpacing','compact','Padding','compact')

% plot the ground truth and estimated source positions
nexttile(1,[num_sources,num_sources])
hold on
scatter(xyb(:,1),xyb(:,2),30*num_sources*xyb(:,3),'k','filled')
scatter(xyb_est(:,1),xyb_est(:,2),30*num_sources*xyb_est(:,3),'r','d','filled')
hold on
legend({'Ground Truth','Estimate'})
axis square; box on; grid on;
xticks(-1:.5:1); yticks(-1:.5:1);
xlim([-1,1]); ylim([-1,1]);
xlabel('$x/\sigma$','interpreter','latex')
ylabel('$y/\sigma$','interpreter','latex')
title('Source Location Estimates','interpreter','latex')

% plot the YKL modes
[X,Y] = meshgrid(sigma*linspace(-5,5,1001));
PSF_modes = 1/sqrt(2*pi*sigma^2)*...
            exp(- (1/(2*sigma)^2)*(...
            (X-permute(xyb_est(:,1),[3,2,1])).^2 +...
            (Y-permute(xyb_est(:,2),[3,2,1])).^2 ));
colormap(turbo)
for t = 1:num_sources
    nexttile(num_sources^2+t)
    mode = sum(permute(A(t,:),[3,1,2]).*PSF_modes,3);
    imagesc([min(X(:)),max(X(:))]/sigma,[min(Y(:)),max(Y(:))]/sigma,abs(mode).^2);
    axis square
    set(gca,'ydir','normal')
    xlabel('$x/\sigma$','interpreter','latex')
    ylabel('$y/\sigma$','interpreter','latex')
    title(sprintf('$|\\pi_{%d}\\rangle$',t),'interpreter','latex')
end

figure
tiledlayout(floor(sqrt(num_sources)),ceil(sqrt(num_sources)),'TileSpacing','compact','Padding','compact')
colormap(turbo)
for t=1:num_sources
    nexttile(t)
    mode = sum(permute(A(t,:),[3,1,2]).*PSF_modes,3);
    %phplot(mode)
    imagesc([min(X(:)),max(X(:))]/sigma,[min(Y(:)),max(Y(:))]/sigma,abs(mode).^2);
    axis square
    set(gca,'ydir','normal')
    xlabel('$x/\sigma$','interpreter','latex')
    ylabel('$y/\sigma$','interpreter','latex')
    title(sprintf('$|\\pi_{%d}\\rangle$',t),'interpreter','latex')

end
end





%% FUNCTIONS 
function Psi = generateIncoherentSourceStateVectors(xyb,sigma)
    % returns the representation of the incoherent source states in the
    % eigenbasis of the density operator.
    % 
    % Psi       : matrix of states. Each column is a unique state.

    sources_xy = xyb(:,1:2);
    priors = xyb(:,3);
    
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
    % Psi       : [K,K] matrix of states (each column is a different state)
    % priors    : [K,1] vector of state priors
    
    % the number of states
    num_states = numel(priors);

    % setup manifold optimization
    problem.M = unitaryfactory(num_states,1);                       % manifold is the space of unitaries (orthnormal projectors)
    problem.cost = @(x) 1 - abs(sum(conj(x).*Psi,1)).^2 * priors;   % cost is the probability of error
    problem.egrad = @(x) - (sum(conj(Psi).*x,1).*priors.').*Psi;  % DO NOT FUCK WITH THIS -- TOOK YOU FOREVER TO DERIVE THE RIGHT EUCLIDEAN GRAD
    
    % solve the optimal measurement
    [YKL,~,~,~] = trustregions(problem);
end

