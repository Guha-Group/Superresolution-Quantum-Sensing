% A script for testing receivers dedicated to multi-source color center
% sensing.
% 
% Author: Nico Deshler
% Date: 4/28/2025

% add utils path
addpath('utils\')

% ground truth parameters
num_sources = 4;
sigma = 1;
xy = sigma*(rand(num_sources,2)-.5);
b = dirchrnd(3*rand(num_sources,1)); 
b = sort(b,'ascend');
xyb = [xy,b];
n_max = 10;

% Photon Allocations
M = 1e6;           % number of calibration photons
N = 1e6;           % number of sensing photons

% simulate SPADE estimation receiver
splitting_ratio = .1;
[n,m] = HGIndices(n_max);
xyb_SS = SimulateReceiver(xyb,M,N,'StaticSPADE',splitting_ratio,[n,m]);

% simulate direct imaging receiver
%[xyb_DI,outdata_DI] = SimulateDirectImaging(xyb,M,N,'DirectImaging',sigma);

function x = dirchrnd(alpha)
    % randomly samples dirichlet random variables with rate parameters
    % given by alpha [K,1]
    y = gamrnd(1,alpha);
    x = y./sum(y,1);
end