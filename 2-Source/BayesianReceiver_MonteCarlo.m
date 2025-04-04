function BayesianReceiver_MonteCarlo(array_id,num_workers,save_dir)
%% COLOR-CENTER SUPERRESOLUTION SENSING MONTE CARLO SIMULATIONS
% Description: 
% Runs a Monte-Carlo analysis for estimating the brightness of
% two sub-diffraction color centers. The analysis steps through different
% pair-wise separations of the two sources with several repeated trials.

% convert the array id input into a numeric index
if ischar(array_id)
    array_id = str2double(array_id);
end

% make a save directory
directory_name = ['data\',save_dir,'\'];
mkdir(directory_name)

% add the utility functions
addpath('utils\')

% setup the monte-carlo survey ranges
sigma = 1;
x0 = 0;
s_range = sigma*(.05:.0125:.5);
k_range = 0:.1:.4;
T = 1; % Trials per monte-carlo sample
M = 5e4;
N = 5e4;
photons_per_adaptation = 1e3;
splitting_ratio = .5; 

% containers for holding results
XSK_0 = zeros(3,numel(s_range),numel(k_range),T);
XSK_EST_DI = XSK_0; 
XSK_EST_SB = XSK_0; 
XSK_EST_AB = XSK_0; 
PHOTONS_AB = zeros(2,numel(s_range),numel(k_range),T);

parpool(num_workers)

% run Monte-Carlo Survey
%for ns = 1:numel(s_range)
    ns = array_id;
    s = s_range(array_id);
    for nk = 1:numel(k_range)
        k = k_range(nk);
        parfor t = 1:T

            % display current iter
            disp([ns,nk,t])

            % random number generation seed should be unique for each
            % worker in parfor loop
            rng(ns+nk+t);

            % construct the parameter vector
            [xsk_0] = [x0,s,k];

            %% Direct Imaging
            [xsk_DI,~] = SimulateReceiver(xsk_0,N,'DirectImaging',M);
            
            %% Static BSPADE
            [xsk_SB,~] = SimulateReceiver(xsk_0,N,'StaticBSPADE',M,splitting_ratio);
            
            %% Adaptive BSPADE
            [xsk_AB,PDF_AB] = SimulateReceiver(xsk_0,N,'AdaptiveBSPADE',M,photons_per_adaptation);
            
            % store results
            XSK_EST_DI(:,array_id,nk,t) = xsk_DI;
            XSK_EST_SB(:,ns,nk,t) = xsk_SB;
            XSK_EST_AB(:,ns,nk,t) = xsk_AB;
            PHOTONS_AB(:,ns,nk,t) = PDF_AB.photons(1:2);

        end
        
        % save 
        save(fullfile(directory_name,['MonteCarlo_2Source_ReceiverSurvey',num2str(array_id),'.mat']),'-regexp', '^(?!ans$).')
    end
%end
end



