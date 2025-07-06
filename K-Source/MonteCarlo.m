function MonteCarlo(array_id, num_workers)
% Runs a Monte-Carlo Survey to determine the performance of different
% estimation algorithms accross a range of scene parameter values
% get array index

% convert array id string to number from slurm
%disp(array_id)
%array_id = str2double(array_id);

% import manopt package
addpath(pwd());

% Recursively add Manopt directories to the Matlab path.
cd('manopt');
addpath(genpath(pwd()));
cd('..');

% add directory paths
addpath('utils/')
addpath('manopt/')

% monte-carlo scene parameter ranges
num_source_range = 3:5;                     % range of number of sub-diffraction sources
min_sep_frac_range = linspace(1/20,1,50);   % range of minimum pair-wise separations to consider
num_trials = 2e3;                           % number of Monte-Carlo Trials per configuration setting
dims = [numel(num_source_range),...
        numel(min_sep_frac_range),...
        num_trials];

% get trial subscripts
[n1,n2,~] = ind2sub(dims,array_id);

% photon allocations
M = 1e6;           % number of calibration photons
N = 1e6;           % number of sensing photons

% receiver/imaging details
sigma = 1;
visualize_flag = 0;
splitting_ratio = 0.1;      % splitting ratio 
n_max = 10;                 % max HG index
[n,m] = HGIndices(n_max);   % HG indices

% number of sources
num_sources = num_source_range(n1);

% minimum separation 
min_sep_frac = min_sep_frac_range(n2); 

% save details
directory = fullfile('MonteCarlo_KSource_v5',[num2str(num_sources),'_sources']);
fname = ['SeparationID_',num2str(n2)];
savefile = fullfile(directory,fname);
mkdir(directory)

% create sub containers
XYB = zeros(num_sources, 3, 1, dims(3));
XYB_SS = zeros(num_sources, 3, 1, dims(3));
XYB_DI = zeros(num_sources, 3, 1, dims(3));

% start parallel loop
%parpool(num_workers);
for n3 = 1:dims(3)  % trials

    % initiate random number gen
    rng(array_id+n3);
    fprintf('Trial %d : %d\n',n3,num_trials)
    
    % generate a scene
    xyb = GenerateRandomConstellationAlt(num_sources,min_sep_frac,ones(num_sources,1),sigma);
    
    % simulate spade receiver
    try
        xyb_SS = SimulateReceiver(xyb,M,N,'StaticSPADE',splitting_ratio,[n,m],sigma,visualize_flag);
    catch
        disp('Error in SPADE receiver')
        xyb_SS = nan(num_sources,3);
    end
    % simulate direct imaging receiver
    xyb_DI = SimulateReceiver(xyb,M,N,'DirectImaging',sigma,visualize_flag);
    
    % store results
    XYB(:,:,1,n3) = xyb;
    XYB_SS(:,:,1,n3) = xyb_SS;
    XYB_DI(:,:,1,n3) = xyb_DI;
end

% save results
save(savefile)
end
