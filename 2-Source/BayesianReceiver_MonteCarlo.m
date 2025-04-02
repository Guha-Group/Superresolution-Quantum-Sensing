%% COLOR-CENTER SUPERRESOLUTION SENSING MONTE CARLO SIMULATIONS
% Description: 
% Runs a Monte-Carlo analysis for estimating the brightness of
% two sub-diffraction color centers. The analysis steps through different
% pair-wise separations of the two sources with several repeated trials.
addpath('utils\')
sigma = 1;
x0 = 0;
s_range = sigma*(.05:.0125:.5);
k_range = 0:.1:.4;
M = 1e4;
N = 5e4;
dm = 1e3;
T = 200; % Trials per monte-carlo sample

% containers for holding MC results
XSK_0 = zeros(3,numel(s_range),numel(k_range),T);
XSK_EST_I = XSK_0;
XSK_EST_II = XSK_0;
XSK_EST_DI = XSK_0;
M_I = zeros(2,numel(s_range),numel(k_range),T);
M_II = M_I;

parpool(14)

% run Monte-Carlo Survey
for ns = 1:numel(s_range)
    s = s_range(ns);
    for nk = 1:numel(k_range)
        k = k_range(nk);
        parfor t = 1:T
            
            % random number generation seed should be unique for each
            % worker in parfor loop
            rng(ns+nk+t);

            % construct the parameter vector
            [xsk_0] = [x0,s,k];

            % simulate the type-I receiver
            [xsk_est_I, PDF_I] = SimulateBayesianReceiver(xsk_0,sigma,M,N,dm,0);

            % simulate the type-II adaptive receiver
            %[xsk_est_II, PDF_II] = SimulateBayesianReceiver(xsk_0,sigma,M,N,dm,1);
            
            % simulate direct imaging receiever
            [xsk_est_DI,~] = SimulateDirectImagingReceiver(xsk_0,sigma,M,N);
            
            % store the results
            XSK_0(:,ns,nk,t) = xsk_0;
            XSK_EST_I(:,ns,nk,t) = xsk_est_I;
            %XSK_EST_II(:,ns,nk,t) = xsk_est_II;
            XSK_EST_DI(:,ns,nk,t) = xsk_est_DI;

            M_I(:,ns,nk,t) = PDF_I.photons(1:2);
            %M_II(:,ns,nk,t) = PDF_II.photons(1:2);

            % display current iter
            disp([ns,nk,t])
           
        end
        
        % save 
        save('MonteCarlo_ReceiverSurvey_TypeI_DI.mat','-regexp', '^(?!ans$).')
    end
end



