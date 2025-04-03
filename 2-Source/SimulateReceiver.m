function [xsk_est, output_data] = SimulateReceiver(xsk_0,N,receiver,varargin)
%% Description:
% This script interfaces a switch block for running a particular receiver
% design on a given input consisting of two color center sources.
%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%
%   xsk_0   :   [1,3] vector of ground truth parameters [x0,s0,k0]
%   N       :   Number of photons allocated to the sensing stage
%   receiver:   string indicating the receiver type
%   varargin:   variable input arguments
%%%%%%%%%%% OUTPUTS %%%%%%%%%%%%%
%   xsk_est     :   [1,3] vector of estimated parameters [x0,s0,k0]
%   output_data :   structure of receiver output metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author(s): Nico Deshler
% date: April 2, 2025

% valid receiver names
validReceivers = {'DirectImaging','StaticBSPADE','AdaptiveBSPADE','ThresholdBSPADE'};

% input validators
validateScene = @(xsk) isnumeric(xsk(1)) & (xsk(2)>0) & (abs(xsk(3))<.5);
validatePositiveInteger = @(x) isnumeric(x) & (rem(x,1)==0) & (x>0);
validatePositiveReal = @(x) isnumeric(x) & (x>0);
validateSplittingRatio = @(x) isnumeric(x) & (0<x) & (x<=1);
validateReceiver = @(x) any(validatestring(x,validReceivers));

% default inputs
default_sigma = 1;

% input parser
p = inputParser;
addRequired(p,'xsk_0',validateScene)
addRequired(p,'N',validatePositiveInteger)
addRequired(p,'receiver',validateReceiver)
addOptional(p,'sigma',default_sigma,validatePositiveReal)

switch receiver
    case 'DirectImaging'
        addRequired(p,'M',validatePositiveInteger)
        parse(p,xsk_0,N,receiver,varargin{:})
        [xsk_est,output_data] = Simulate_DirectImaging_Receiver(xsk_0,N,p.Results.M,p.Results.sigma);

    case 'StaticBSPADE'
        addRequired(p,'M',validatePositiveInteger)
        addRequired(p,'splitting_ratio',validateSplittingRatio)
        parse(p,xsk_0,N,receiver,varargin{:})
        [xsk_est,output_data] = Simulate_StaticBSPADE_Receiver(xsk_0,N,p.Results.M,p.Results.sigma,p.Results.splitting_ratio);
    
    case 'AdaptiveBSPADE'
        addRequired(p,'M',validatePositiveInteger)
        addRequired(p,'dm',validatePositiveInteger)
        parse(p,xsk_0,N,receiver,varargin{:})
        [xsk_est,output_data] = Simulate_AdaptiveBSPADE_Receiver(xsk_0,N,p.Results.M,p.Results.sigma,p.Results.dm);
    
    case 'ThresholdBSPADE'
        addRequired(p,'dm',validatePositiveInteger)
        addRequired(p,'target_std_s',validatePositiveReal)
        parse(p,xsk_0,N,receiver,varargin{:})
        [xsk_est,output_data] = Simulate_ThresholdBSPADE_Receiver(xsk_0,N,p.Results.sigma,p.Results.dm,p.Results.target_std_s);
end
end