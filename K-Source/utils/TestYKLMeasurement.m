% This script explores the use of manopt for determining the optimal
% Helstrom measurement to apply for discriminating quantum states. All
% examples are for particular point source configurations.
%
% Nico Deshler

%addpath('.\manopt\*')


%% K POINT SOURCES

% preliminaries
sigma = 1;                                  % Rayleigh Limit (for gaussian PSF)
num_sources = 3;                            % number of sources in scene
xy = sigma*randn([num_sources,2]);          % source cartesian coordinates
priors = permute([1/3,1/3,1/3]',[3,2,1]);   % brightness (priors) for each source state
[states,Psi] = generateIncoherentSourceStateVectors(xy,priors,sigma);
rho = sum(priors.*states,3);                % density operator

% evaluate the SLD measurements
L = getSLD(states,priors);

% get the mean SLD
L0 = sum(priors(:,:,1:2).*L,3);

% determine expected commutator of the SLDs (sufficient condition for
% simultaneous QFIM saturability)
avg_commutators = zeros((num_sources-1)*[1,1]);
for j = 1:(num_sources-1)
    for k = 1:(num_sources-1)
        avg_commutators(j,k) = trace(rho*L(:,:,j)*L(:,:,k));
    end
end

% determine the realness of the SLDs commutators 
% If all are real then the Helstrom Cramer Rao
% bound is equal to the Quantum Cramer Rao bound.
HCRB_equals_QCRB = all(imag(avg_commutators(:))==0);

% evaluate the minimum error probability measurement (Helstrom measurement)
M = getYKL(states,priors);
W_M = POVM2projectors(M); % vectorize the POVM projectors into column vectors
M_alt = (YKLMeasurement(Psi'*Psi,diag(squeeze(priors)))/(Psi))';
W_L = M_alt;

% vectorize the SLD into column vectors
%W_L = SLD2projectors(L);
W_L0 = SLD2projectors(L0);



%plot2DProjectors(W_M,W_L,Psi,priors)
plot3DProjectors(W_M,W_L,W_L0,Psi,priors)



%% FUNCTIONS 
function [states,Psi] = generateIncoherentSourceStateVectors(sources_xy,priors,sigma)
    % returns the representation of the incoherent source states in the
    % eigenbasis of the density operator.
    % 
    % states    : pure states as projectors (stack of dyad matrices)
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

    % construct the states as density operators
    states = permute(Psi,[1,3,2]).*permute(conj(Psi),[3,1,2]);

end


function M = getYKL(states,priors)
    % Compute the Helstrom/YKL measurement for single-shot state discrimination via
    % optimization of the error probability.

    % the number of states
    num_states = numel(priors);

    % setup manifold optimization
    manifold = sympositivedefinitesimplexcomplexfactory(num_states,num_states);
    problem.M = manifold;

    % cost function is the negative of the probability of correct detection
    problem.cost = @(x)  - sum(priors.* real(sum(pagetranspose(states).*x,[1,2])),3); % equivalent to sum_i p_i trace(x_i rho_i))
    problem.egrad = @(x) - priors.*pagetranspose(states);

    % solve the optimization
    % solve the optimal measurement
    [M,~,~,~] = trustregions(problem);

    % check that the YKL conditions are satisfied.
    C1 = zeros(num_states*[1,1,1,1]);
    % condition (1): All elements in C1 should be zero
    for j = 1:num_states
        for k = 1:num_states           
            C1(:,:,j,k) = M(:,:,j)*(priors(j)*states(:,:,j) - priors(k)*states(:,:,k))*M(:,:,k);
        end

    end 
    % condition (2): All matrices in C2 should be PSD
    C2 = sum(priors.*pagemtimes(states,M),3) - priors.*states;
end


function L = getSLD(states,priors)
    % solves the Lyapunov equation under assumption that the parameters of
    % interest are the priors
    num_states = numel(priors);
    L = zeros(size(states)-[0,0,1]);
    rho = sum(priors.*states,3);
    for k = 1:(num_states-1)
        L(:,:,k) = lyap(rho,-2*(states(:,:,k)-states(:,:,end)));
  %      L(:,:,k) = lyap(rho,-2*states(:,:,k));
    end
end

function T = getPGM(states,priors)
    % returns the "pretty good measurement"
    num_states = numel(priors);
    T = zeros(size(states));
    

end

function S = getSIC_POVM(n)
    % returns the symmetric information complete (SIC) POVM for the hilbert
    % space of dimension n. This is for quantum state tomography.
end

function W = POVM2projectors(M)
    % converts a stack of POVM operators to projectors where each column of
    % W is a pure-state vector.
    W = pagemtimes(M,ones([size(M,1),1]));
    W = W./ sqrt(sum(abs(W).^2,1));
    W = permute(W,[1,3,2]);
end

function W = SLD2projectors(L)
% converts a stack of SLDs into a stack of projectors where each column of
% W is a pure-state vector.
W = zeros(size(L));
for k = 1:size(L,3)
    [U_L,D_L] = eig(L(:,:,k)); % get the eigenvectors of the SLD operator
    W(:,:,k) = U_L;
end

end

function plot2DProjectors(W_L,W_M,psi,priors)
    % plots the numerically-optimal 2D projection vectors for binary state
    % discrimination with arbitrary priors and inner products. This script
    % shows both the projectors for the Helstrom measurement and the SLD
    % measurement associated with estimating the prior probability.
    
    % get the dimensionality of the space
    hilbert_dim = size(W_M,2);
    assert(hilbert_dim==2);

    % some coordinates for plotting a reference circle
    [xc,yc] = pol2cart(linspace(0,2*pi,1001),ones(1,1001)); % reference unit circle
    
    % negate the projector such that the Helstrom measurement lies in the 
    % positive half of the cartesian plane and the SLD measurement
    % lies in the bottom half (more for aesthetics)
    W_M = +W_M.*sign(W_M(2,:));
    W_L = -W_L.*sign(W_L(2,:));

    % plot the state and projection vectors for optimal measurement
    figure
    hold on
    % plot states
    quiver(0,0,psi(1,1),psi(2,1),0,'k','LineWidth',1.5)
    quiver(0,0,psi(1,2),psi(2,2),0,'k','LineWidth',1.5)
    % plot helstrom measurement
    quiver(0,0,W_M(1,1),W_M(2,1),0,'r','LineWidth',1.5)
    quiver(0,0,W_M(1,2),W_M(2,2),0,'r','LineWidth',1.5)
    % plot SLD measurement
    quiver(0,0,W_L(1,1),W_L(2,1),0,'b','LineWidth',1.5)
    quiver(0,0,W_L(1,2),W_L(2,2),0,'b','LineWidth',1.5)
    % plot reference lines
    plot(xc,yc,'--k');
    xline(0,'k','LineWidth',.5);
    yline(0,'k','LineWidth',.5);
    hold off
    axis square
    box on
    xticks((-1:.5:1));
    yticks((-1:.5:1));
    xlim([-1,1])
    ylim([-1,1])
    xlabel('$\langle e_{-} |\psi \rangle$','interpreter','latex')
    ylabel('$\langle e_{+} |\psi \rangle$','interpreter','latex')
    legend({'$|\psi_1\rangle$','$|\psi_2\rangle$',...
            '$|\pi_1\rangle$ (Helstrom)','$|\pi_2\rangle$ (Helstrom)',...
            '$|\pi_1\rangle$ (SLD)','$|\psi_2\rangle$ (SLD)'},'interpreter','latex')
    title(sprintf('Prior $p_1 : %.2f$',priors(1)),'interpreter','latex')
end


function plot3DProjectors(W_M,W_L,W_L0,psi,priors)
    % plots the numerically-optimal 3D projection vectors for ternary state
    % discrimination with arbitrary priors and inner product. This function
    % shows both the projectors for the Helstrom measurement and the SLD
    % measurement associated with estimating the prior probability.
    
    % get the dimensionality of the space
    hilbert_dim = size(W_M,2);
    assert(hilbert_dim==3);

    % negate the projectors such that the Helstrom measurement lives in the
    % upper hemisphere of the 3D Ball
    psi = +psi.*sign(psi(3,:));
    W_M = +W_M.*sign(W_M(3,:));
    W_L = +W_L.*sign(W_L(3,:,:));
    W_L0 = +W_L0.*sign(W_L0(3,:));

    % plot the state and projection vectors for optimal measurement
    figure
    hold on
    % plot states
    for n = 1:hilbert_dim
        quiver3(0,0,0,psi(1,n),psi(2,n),psi(3,n),0,'k','LineWidth',1.5)
    end

    % plot helstrom measurement
    for n=1:hilbert_dim
        quiver3(0,0,0,W_M(1,n),W_M(2,n),W_M(3,n),0,'r','LineWidth',1.5)
    end
    
    % plot combined SLD measurement projectors
    for n=1:hilbert_dim
        quiver3(0,0,0,W_L0(1,n),W_L0(2,n),W_L0(3,n),0,'b','LineWidth',1.5)
    end

    % plot SLD measurements
    sld_colors = turbo(size(W_L,3));
    for m = 1:size(W_L,3)
        for n=1:hilbert_dim
            quiver3(0,0,0,W_L(1,n,m),W_L(2,n,m),W_L(3,n,m),0,'--','Color',sld_colors(m,:),'LineWidth',1.5)
        end
    end
    
    legend({'Pure States','','','Helstrom Projectors','','','Combined SLD'})
    
    hold off
    axis square
    box on
    grid on
    xlim([-1,1]);    ylim([-1,1]);    zlim([-1,1]);

    % plot the probability of error for each measurement
    figure
    P_err_M =  1-squeeze(priors)'*diag(abs(W_M'*psi).^2);
    P_err_L0 = 1-squeeze(priors)'*diag(abs(W_L0'*psi).^2);
    P_err_M_alt = 1-squeeze(priors)'*diag(abs(W_L'*psi).^2);
    
    bar([P_err_L0,P_err_M,P_err_M_alt]);
    xticklabels({'Combined SLD','Helstrom (YKL)','YKL (alt)'})
    ylabel('$P_{e}$','interpreter','latex')
    xlabel('Measurements','interpreter','latex')

end


%{
%% TWO POINT SOURCES
num_states=2;
s = 1;
sigma = 1;
g = exp(-(s/sigma)^2); % inner product between two point sources

% 2D vector for emitter states in |e_+> |e_-> representation
v1 =  [sqrt(1+g),sqrt(1-g)]'/sqrt(2);
v2 = [sqrt(1+g),-sqrt(1-g)]'/sqrt(2);
states(:,:,1) = v1*v1';
states(:,:,2) = v2*v2';
p = .25;
priors = permute([p,1-p]',[3,2,1]);

manifold = sympositivedefinitesimplexcomplexfactory(num_states,num_states);
problem.M = manifold;


% cost function is the negative of the probability of correct detection
problem.cost = @(x)  - sum(priors.* real(sum(pagetranspose(states).*x,[1,2])),3); % equivalent to sum_i p_i trace(x_i rho_i))
problem.egrad = @(x) - priors.*pagetranspose(states);

% numerically check gradient consistency
%checkgradient(problem);
%pause;

% solve the optimal measurement
[x,xcost,info,options] = trustregions(problem);
%}
