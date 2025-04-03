%% Test Receivers
% A simple test script to ensure all of the receivers are working without
% throwing any errors.
%
% Author(s): Nico

% add utils path
addpath('utils\')

% ground truth parameters
xsk_0 = [0,.01,-.2]; % [midpoint,half-sep,brightness bias]
M = 1e5;           % number of calibration photons
N = 5e4;           % number of sensing photons

%% Direct Imaging
[xsk_DI,PDF_DI] = SimulateReceiver(xsk_0,N,'DirectImaging',M);

%% Static BSPADE
splitting_ratio = 0.5;
[xsk_SB,PDF_SB] = SimulateReceiver(xsk_0,N,'StaticBSPADE',M,splitting_ratio);

%% Adaptive BSPADE
photons_per_adaptation = 1e3;
[xsk_AB,PDF_AB] = SimulateReceiver(xsk_0,N,'AdaptiveBSPADE',M,photons_per_adaptation);



%% Threshold BSPADE
% photons_per_adapatation = 1e3; target_std_s=1e-3;
%[xsk_TB,PDF_TB] = SimulateReceiver(xsk_0,N,'ThresholdBSPADE',photons_per_adapatation,target_std_s);





