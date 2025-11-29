function [XPR,searchRad] = sim_and_apply_DETRIM(D,locPrec,searchTime,varargin)
% SIM_AND_APPLY_DETRIM simulates a immobile/mobile particle and executes the DETRIM
% clustering algorithm in a single-pass mode to calculate the fraction of
% points incorrectly classified as immobile.
%
% This function is primarily used to determine the False Positive Rate (FPR)
% for mobile particles, where a high FPR indicates DETRIM is confusing
% mobile diffusion events with true immobile events.
%
%   INPUTS:
%   D: scalar; Diffusion coefficient [px^2/frame].
%   locPrec: scalar; Localization precision [px].
%   searchTime: scalar; Time search window [frames].
%
%   OUTPUTS:
%   XPR: scalar; The True/False Positive Rate.
%
%   Copyright (c) 2025 Christian Paolo Richter
%   University of Osnabrueck

ip = inputParser;
addRequired(ip,'D')
addRequired(ip,'locPrec')
addRequired(ip,'searchTime')
addParamValue(ip,'T',100000) % Total number of frames (simulation duration)
addParamValue(ip,'P',0.99) % Expected True Positive Rate for calculating the search radius
parse(ip,D,locPrec,searchTime,varargin{:});

% 1. Simulate a 2D isotropic Brownian trajectory (Mobile, D > 0)
[x,y,t] = ...
    sim_2D_isotropic_brownian_trajectory(...
    1,ip.Results.T,D,locPrec);

% 2. Calculate the necessary search radius (rSearch)
% The formula ensures P*100% probability of re-detection for an immobile particle
searchRad = 2*locPrec*sqrt(-log(1-ip.Results.P)); %[px]

% 3. Run DETRIM with single-pass clustering
clusterID = ...
    DETRIM(...
    x,y,t,...
    searchRad,...
    searchTime);

% 4. Evaluation: Calculate total points per cluster (N)
N = accumarray(clusterID,1);

% Fraction unclassified points (N(1) being the noise cluster ID 1)
pctUnclass = N(1)/ip.Results.T;

%XPR equals TPR in case of immobile and FPR in case of mobile particles simulated
XPR = 1-pctUnclass;
end %fun