% A script for testing receivers dedicated to multi-source color center
% sensing.
% 
% Author: Nico Deshler
% Date: 4/28/2025

% add utils path
addpath('utils\')

% ground truth parameters
sigma = 1;                          
num_sources = 4;
min_sep_frac = sigma/8;
contrast = .1;
visualize_flag = 1;

% generate scene
xyb = GenerateRandomConstellation(num_sources,min_sep_frac,contrast,sigma);

% Photon Allocations
M = 1e5;           % number of calibration photons
N = 1e5;           % number of sensing photons

% simulate SPADE estimation receiver
splitting_ratio = .1;
n_max = 10;
[n,m] = HGIndices(n_max);
xyb_SS = SimulateReceiver(xyb,M,N,'StaticSPADE',splitting_ratio,[n,m],sigma,visualize_flag);

% simulate direct imaging receiver
xyb_DI = SimulateReceiver(xyb,M,N,'DirectImaging',sigma,visualize_flag);
