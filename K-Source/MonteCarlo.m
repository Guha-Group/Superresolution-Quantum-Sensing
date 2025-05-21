% Runs a Monte-Carlo Survey to determine the performance of different
% estimation algorithms accross a range of scene parameter values

% add utils path
addpath('utils\')

% save filename
fname = 'MonteCarlo_KSource_SPADE_v_DirectImaging.mat';

% monte-carlo scene parameter ranges
num_source_range = 3:5;                         % range of number of sub-diffraction sources
min_sep_frac_range = linspace(sigma/20,1,50);   % range of minimum pair-wise separations to consider
contrast_range = .5:.1:1;       % range of contrasts to explore
num_trials = 200;               % number of Monte-Carlo Trials per configuration setting
dims = [numel(num_source_range),numel(min_sep_frac_range),numel(contrast_range),num_trials];

% photon allocations
M = 1e5;           % number of calibration photons
N = 1e5;           % number of sensing photons

% receiver/imaging details
sigma = 1;
visualize_flag = 0;
switching_ratio = 0.1;      % switching ratio 
max_n = 10;                 % max HG index
[n,m] = HGIndices(n_max);   % HG indices

% make data containers
XYB = cell(dims(1),1);
XYB_SS = cell(dims(1),1);
XYB_DI = cell(dims(1),1);

for n1 = 1:dims(1)    % number of sources
    num_sources = num_source_range(n1);

    % create sub containers
    XYB_n1 = zeros(num_sources, 3, dims(2),dims(3),dims(4));
    XYB_SS_n1 = zeros(num_sources, 3, dims(2),dims(3),dims(4));
    XYB_DI_n1 = zeros(num_sources, 3, dims(2),dims(3),dims(4));

    for n2 = 1:dims(2)  % contrast range
        contrast = contrast_range(n2);
        for n3 = 1:dims(3)  % separation range
            min_sep_frac = min_sep_frac_range(n3); 
            for n4 = 1:dims(4)  % trials

                % generate a scene
                xyb = GenerateRandomConstellation(num_sources,min_sep_frac,contrast,sigma);
                
                % simulate spade receiver
                xyb_SS = SimulateReceiver(xyb,M,N,'StaticSPADE',splitting_ratio,[n,m],sigma,visualize_flag);
                
                % simulate direct imaging receiver
                xyb_DI = SimulateReceiver(xyb,M,N,'DirectImaging',sigma,visualize_flag);

                % store results
                XYB_n1(:,:,n2,n3,n4) = xyb;
                XYB_SS_n1(:,:,n2,n3,n4) = xyb_SS;
                XYB_DI_n1(:,:,n2,n3,n4) = xyb_DI;
            end
            % save intermediate results
            save(fname)
        end
    end
    % store complete data arrays
    XYB{n1} = XYB_n1;
    XYB_SS{n1} = XYB_SS_n1;
    XYB_DI{n1} = XYB_DI_n1;

    % save general results
    save(fname)
end
