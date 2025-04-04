function BayesianReceiver_MonteCarlo(array_id,num_workers,save_dir)
%% COLOR-CENTER SUPERRESOLUTION SENSING MONTE CARLO SIMULATIONS
% Description: 
% Runs a Monte-Carlo analysis for estimating the brightness of
% two sub-diffraction color centers. The analysis steps through different
% pair-wise separations of the two sources with several repeated trials.

% convert string inputs into numerics
if ischar(array_id)
    array_id = str2double(array_id);
end
if ischar(num_workers)
    num_workers = str2double(num_workers);
end

% make a save directory
directory_name = ['data/',save_dir,'/'];
mkdir(directory_name)

% add the utility functions
addpath('utils/')

% setup the monte-carlo survey ranges
sigma = 1;
x0 = 0;
k_range = 0:.1:.4;
s_range = sigma*linspace(.01,.5,50);
T = 500; % Trials per monte-carlo sample
M = 1e4;
N = 1e4;
photons_per_adaptation = 1e3;
splitting_ratio = .5; 


% collect all simulation variables and meta-data into a structure
sim_struct.directory_name = directory_name;
sim_struct.array_id = array_id;
sim_struct.sigma = sigma;
sim_struct.x0 = x0;
sim_struct.k_range = k_range;
sim_struct.s_range = s_range;
sim_struct.T = T;
sim_struct.M = M;
sim_struct.N = N;
sim_struct.photons_per_adaptation = photons_per_adaptation;
sim_struct.splitting_ratio = splitting_ratio;


% containers for holding results
XSK_0 = zeros(3,numel(s_range),numel(k_range),T);
XSK_EST_DI = XSK_0; 
XSK_EST_SB = XSK_0; 
XSK_EST_AB = XSK_0; 
PHOTONS_AB = zeros(2,numel(s_range),numel(k_range),T);

parpool(num_workers)

% run Monte-Carlo Survey
nk = array_id;
k = k_range(nk);
for ns = 1:numel(s_range)
    s = s_range(ns);
    parfor t = 1:T
    
        % display current iter
        disp([nk,ns,t])
    
        % random number generation seed should be unique for each
        % worker in parfor loop
        rng(nk+ns+t);
    
        % construct the parameter vector
        xsk_0 = [x0,s,k];
    
        %% Direct Imaging
        [xsk_DI,~] = SimulateReceiver(xsk_0,N,'DirectImaging',M);
        
        %% Static BSPADE
        [xsk_SB,~] = SimulateReceiver(xsk_0,N,'StaticBSPADE',M,splitting_ratio);
        
        %% Adaptive BSPADE
        [xsk_AB,PDF_AB] = SimulateReceiver(xsk_0,N,'AdaptiveBSPADE',M,photons_per_adaptation);
        
        % store results
        XSK_0(:,nk,ns,t) = xsk_0;
        XSK_EST_DI(:,nk,ns,t) = xsk_DI;
        XSK_EST_SB(:,nk,ns,t) = xsk_SB;
        XSK_EST_AB(:,nk,ns,t) = xsk_AB;
        PHOTONS_AB(:,nk,ns,t) = PDF_AB.photons(1:2);
    
    end
    
    % save 
    save(fullfile(directory_name,['MonteCarlo_2Source_ReceiverSurvey',num2str(array_id),'.mat']),...
        'sim_struct','XSK_0','XSK_EST_DI','XSK_EST_SB','XSK_EST_AB','PHOTONS_AB')
end
%end
end



