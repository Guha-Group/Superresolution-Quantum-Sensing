function [xsk_est,outdata] = Simulate_DirectImaging_Receiver(xsk_0,N,M,sigma)
% Description: Simulates a two-stage receiver consisting only of direct
% imaging measurements.

%%%%%%%%%%% INPUTS %%%%%%%%%%%%%
%   xsk_0   :   [1,3] vector of ground truth parameters [x0,s0,k0]
%   M       :   Number of photons allocated to Stage 1. If the
%               adaptive flag is false, all M photons are allocated to
%               the BSPADE measurement and we integrate as many photons
%               as needed in the direct imaging phase of the calibration 
%               to guarantee good performance with the Binary SPADE 
%               stage. If the adaptive flag is true, then the M photons
%               are adaptively allocated between Direct Imaging and BSPADE
%   N       :   Number of photons allocated to Stage 2.
%
%%%%%%%%%%% OUTPUTS %%%%%%%%%%%%
%   xsk_est :   [1,3] vector of parameter estimates [x0_mle,s0_mmmse,k0_mmse]
%   PDF     :   A structure with all priors and marginalized posterior
%               distributions of the estimators and their domains.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author(s):    Nico Deshler, ndeshler@arizona.edu
% institution:  University of Arizona, Wyant College of Optical Sciences
% date:         April 1, 2025

% just return a static BSPADE receiver with all photons allocated to the
% direct imaging stage
splitting_ratio = 1;
[xsk_est, outdata] = Simulate_StaticBSPADE_Receiver(xsk_0,N,M,sigma,splitting_ratio);

end