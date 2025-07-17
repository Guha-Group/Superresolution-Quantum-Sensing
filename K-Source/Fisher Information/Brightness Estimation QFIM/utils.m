classdef utils
    methods (Static)

        function [simplex_samples,num_simplex_samples] = sample_3_simplex(num_samples_per_dim)
            % Generates all combinations of possible weight vectors on the 3-simplex 
            % via a 'stick-breaking' method. The output is an array is a 3xM matrix
            % where each each column is a different sample in the probability simplex.
            % (i.e. each column sums to 1 and is greater than 0 - edges are excluded)
        
            q = linspace(0,1,num_samples_per_dim);
            num_simplex_samples = (numel(q)+1)*numel(q)/2;
            simplex_samples = zeros(3,num_simplex_samples);
            k=1;
            for i=1:numel(q)
                for j=i:numel(q)
                    p1 = q(i); p2 = q(j)-q(i); p3 = 1-p1-p2;
                    simplex_samples(:,k) = [p1,p2,p3]';
                    k = k+1;
                end
            end
        
            % prune samples that lie on boundary of simplex
            simplex_samples = simplex_samples(:,prod(simplex_samples,1)>0);
            num_simplex_samples = size(simplex_samples,2);
        end
        
        function psi = scene_states(xyb)
            % Computes the representation of the incoherent point-source states 
            % in the eigenbasis of the density operator.
            % 
            % INPUT
            % xyb       : [Kx3] real matrix. Each row is (x_i,y_i,b_i) for
            %                   a scene consisting of K point sources.
            %
            % OUTPUT
            % psi       : [KxK] complex matrix. Each column is a unique state.
            
            xy = xyb(:,1:2);
            b = xyb(:,3);
            
            % compute the pairwise difference vectors between all source positions
            deltas_xy = permute(xy,[1,3,2])-permute(xy,[3,1,2]);
            
            % make diagonal matrix of prior probabilities
            P = diag(b(:));
        
            % Gram Matrix (assuming gaussian PSF)
            G = exp(-sum(deltas_xy.^2,3)/8);
          
            % Get spectral decomposition of Gram matrix
            [U,D] = eig(G);
            
            % Determine the S matrix
            A = sqrt(P)*U*sqrt(D);
            S = A'*A;
            
            % Get spectral decomposition of S matrix (the eigenvalues of S are the same
            % as the eigenvalues of rho)
            [V_dagger, rho_evals] = eig(S);
            V = V_dagger';
            
            % Compute the representation matrix of psi in the 
            % eigenbasis of rho
            psi = V * sqrt(D) * U';
        end
        
        
        function [ykl, p_e] = ykl_states(states, priors)
            % computes the Yuen-Kennedy-Lax projectors that optimally distinguish a
            % set of 'states' with 'priors'.
            %
            % INPUTS:
            % states    :   [K x K] complex matrix where each column is a unique state
            % priors    :   [K x 1] prior probability vector where all entries sum
            %                       to 1
            % OUTPUTS:
            % ykl       :   [K x K] complex matrix of YKL projectors where each
            %                       column is a unique state. 'ykl' is unitary.
            % p_e       :   [1 x 1] the state classification error probability 
            %                        achieved under the ykl measurement
            
            % dimensionality of state space
            K = size(states,1);
        
            % manifold of unitary matrices
            manifold = unitaryfactory(K, 1);
            
            % cost function - probability of error
            P_err = @(U) 1 - trace(diag(priors).*abs(U'*states).^2);
            cost = @(U) P_err(U);
            
            problem.M = manifold;
            problem.cost = cost;
            options.verbosity = 0;
            warning('off', 'manopt:getGradient:approx');
            warning('off', 'manopt:getHessian:approx');
            
            % more magic solver
            [ykl, p_e, ~] = trustregions(problem,[],options);
        end

        function Q = brightness_params_QFIM(states, priors)
            % Returns the QFI matrix for the priors (treated as parameters) 
            % in a density operator consisting of K dyads. That is,
            % rho = 
            % The last brightness parameter b_K = 1 - (b_1 + ... + b_[K-1])
            %
            % INPUTS:
            % states    :   [K x K] complex matrix where each column is a unique state
            % priors    :   [K x 1] prior probability vector where all entries sum
            %                       to 1
            %
            % OUTPUTS:
            % Q         :   [K-1 x K-1] state matrix

            % dimensionality of space
            K = size(states, 1);

            % density operator
            P = diag(priors);
            rho = states*P*states';

            % instantiate QFI matrix
            Q = zeros(K-1);
            for i = 1:(K-1)
                
                % partial w.r.t brightnesss b_i
                drho_di = states(:,i)*states(:,i)' - states(:,end)*states(:,end)';
                
                % SLD for b_i
                L_i = lyap(rho,-2*drho_di);
                
                for j = 1:(K - 1)
                    % partial w.r.t brightnesss b_j
                    drho_dj = states(:,j)*states(:,j)' - states(:,end)*states(:,end)';

                    % SLD for b_j
                    L_j = lyap(rho,-2*drho_dj);
                    
                    % QFI matrix element
                    Q(i,j) =  trace(rho*(L_i *L_j + L_j * L_i)) / 2;
                end
            end
        end

        function [Q,C] = tangent_brightness_QFIM(states,priors)
            % returns the QFIM and the expected commutation relation of the
            % SLDs

            % number of states
            K = size(states,1);

            % tangent space matrix for the brightness perturbation space
            V = eye(K) - ones(K,K)/K;

            % density operator
            P = diag(priors);
            rho = states*P*states';

            % make vector of states
            rho_vec = permute(states,[1,3,2]).*permute(conj(states),[3,1,2]);

            % get perturbation states
            nu = permute(pagemtimes(V,permute(rho_vec,[3,1,2])),[2,3,1]);
            
            % get SLDs
            L = zeros(K,K,K);
            for k = 1:K
                L(:,:,k) = lyap(rho,-2*nu(:,:,k));
            end
            
            % define the inner product and anti-product operations
            ip = @(A,B) trace(rho*(A*B + B*A))/2;
            ap = @(A,B) trace(rho*(A*B - B*A))/2;

            % get commutation matrix
            C = zeros(K,K);
            for i=1:K
                for j=1:K
                    C(i,j) = ap(L(:,:,i),L(:,:,j));
                end
            end

            % get QFIM
            Q = zeros(K,K);
            for i = 1:K
                for j = 1:K
                    Q(i,j) = ip(L(:,:,i),L(:,:,j));
                end
            end

        end
        

        function [cfim_numeric,cfim_analytic] = pgm_cfim(states,priors)
           
            % number of states
            K = size(states,1);

            % tangent space matrix for the brightness perturbation space
            V = eye(K) - ones(K,K)/K;
            
            % get ground truh tangent params
            theta_0 = pinv(V)*(priors-1/K);

            % make vector of states
            rho_vec = permute(states,[1,3,2]).*permute(conj(states),[3,1,2]);
            
            % compute the pgm
            rho = @(theta) mean(rho_vec,3) + sum(permute(V*theta,[3,2,1]).*rho_vec,3);
            rho_insq = @(theta) inv(sqrtm(rho(theta)));
            pgm = @(theta) pagemtimes(rho_insq(theta),...
                  pagemtimes(...
                  permute(1/K + V*theta,[3,2,1]).*rho_vec,rho_insq(theta)));
            

            % reference density operator (for debugging)
            rho_ref = sum(permute(priors,[3,2,1]).*rho_vec,3);
            rho_0  = rho(theta_0);


            % compute measurement probabilities
            p = @(theta) sum(pagemtimes(rho(theta),pgm(theta)).*eye(K),[1,2]);
            
            % calculate the CFIM via numerical differentiation
            cfim_numeric = zeros(K);
            I = eye(K);
            eps = 1e-10;
            p0 = p(theta_0); % base probability
            for i = 1:K
                for j = 1:K
                    dtheta_i = eps*I(:,i);
                    dtheta_j = eps*I(:,j);
                    dp_i = (p(theta_0 + dtheta_i) - p0)/eps;
                    dp_j = (p(theta_0 + dtheta_j) - p0)/eps;
                    cfim_numeric(i,j) = sum(dp_i.*dp_j./p0,3);
                end
            end

            % get analytic cfim
            cfim_analytic = V*diag(1./priors)*V.';

        end

        function xy_t = bezier(xy_c,num_pts)
            % returns num_pts coords along a bezier curve defined by the control
            % points xy_c. A beautiful recursive definition.
            
            t = linspace(0,1,num_pts);
            xy_t = zeros([num_pts,size(xy_c,[2,3])]);
            for n = 1:num_pts
                xy_t(n,:) = bezier_point(xy_c,t(n));
            end
            
            function xy_t = bezier_point(xy, t)
                if size(xy,1) == 1
                    xy_t = xy;
                else 
                    xy_t = (1-t)*bezier_point(xy(1:end-1,:),t)+t*(bezier_point(xy(2:end,:),t));
                end
            end
        end


    end
end