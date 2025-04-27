%% Test Receivers
% A simple test script to ensure all of the receivers are working properly.
%
% Author(s): Nico
% Date: April 26, 2025

% add utils path
addpath('utils\')

% ground truth parameters
xsk_0 = [0,.01,.4]; % [midpoint,half-sep,brightness bias]
M = 1e5;           % number of calibration photons
N = 1e5;           % number of sensing photons

%% Direct Imaging
%[xsk_DI,outdata_DI] = SimulateReceiver(xsk_0,N,'DirectImaging',M);
%PlotReceiver(outdata_DI,'DirectImaging');

%% Static BSPADE
%splitting_ratio = .5;
%[xsk_SB,outdata_SB] = SimulateReceiver(xsk_0,N,'StaticBSPADE',M,splitting_ratio);
%PlotReceiver(outdata_SB,'StaticBSPADE');

%% Adaptive BSPADE
photons_per_adaptation = 1e3;
[xsk_AB,PDF_AB] = SimulateReceiver(xsk_0,N,'AdaptiveBSPADE',M,photons_per_adaptation);
PlotReceiver(PDF_AB,'AdaptiveBSPADE')

%% Threshold BSPADE
% photons_per_adapatation = 1e3; target_std_s=1e-3;
%[xsk_TB,PDF_TB] = SimulateReceiver(xsk_0,N,'ThresholdBSPADE',photons_per_adapatation,target_std_s);
%PlotReceiver(PDF_tb,'ThresholdBSPADE')