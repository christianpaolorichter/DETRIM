function [clusterID,pntType] = DETRIM(x,y,t,searchRad,searchTime,varargin)
% DETRIM (DEtection of TRansient IMmobilization) is a robust, hierarchical
% DBSCAN-based algorithm designed to detect both long-lived and transient
% immobilization events in Single-Molecule Localization Microscopy (SMLM) data.
%
% This function adapts the DBSCAN principle of density-reachability over a
% range of time windows (hierarchical search) to group subsequent observations
% of the same immobile single-molecule. It uses forward and reverse clustering
% followed by fusion to mitigate temporal edge effects and ensure complete
% capture of fragmented or blinking events.
%
%   INPUTS:
%   t: vector; timepoint (image frame) at which the molecule has been localized
%   x: vector; x-position of the molecule [px]
%   y: vector; y-position of the molecule [px]
%   searchRad: scalar; defines the max. distance [px] between observations to be considered as potential links
%   searchTime: defines the temporal search window  and the minimum duration (in frames) required for an immobile event.
%               If scalar: Defines a single temporal search window (in frames) AND the minimum number of sequential frames required to define a transient immobile event.
%               If vector: Defines multiple temporal search windows (in frames), in ascending order, used for the hierarchical search. 
%               The first element also defines the minimum number of sequential frames required to define a transient immobile event.
% 
% 
%   OPTIONAL PARAMETERS:
%   xFactor: scalar; (DEFAULT: 1.0). The Point Density Correction Factor. Must be between 0 and 1.
%            Used to account for fluorophore blinking/missing data by defining the
%            minimum required number of points for core classification as: (searchTime * xFactor).
%   verbose: scalar; (DEFAULT: 0). If 1, calls DETRIM_verbose to plot the resulting clusters.
%
%   OUTPUTS:
%   clusterID:	vector; Unique identifier (integer) for each detected immobile cluster. ID=1 is reserved for mobile/noise points.
%   pntType	vector; Classification of each point: 1=Core Point, 0=Border Point, -1=Noise/Mobile Point.
%
%   Copyright (c) 2025 Christian Paolo Richter
%   University of Osnabrueck
%

%%
ip = inputParser;
ip.KeepUnmatched = true;
addRequired(ip,'x')
addRequired(ip,'y')
addRequired(ip,'t')
addRequired(ip,'searchRad', @(x)isscalar(x) && x > 0)
addRequired(ip,'searchTime', @(x)(isscalar(x) || (isvector(x) && issorted(x))) && all(x > 0))
addParamValue(ip,'xFactor',1, @(x)isscalar(x) && x > 0 && x <= 1) %#ok
addParamValue(ip,'verbose',0, @(x)isscalar(x)) %#ok
parse(ip,x,y,t,searchRad,searchTime,varargin{:});

xFactor = ip.Results.xFactor;
verbose = ip.Results.verbose;

%% CLUSTERING
%% iterative segmentation of the clusters
take = 1:numel(t); %we start with all points in the game
clusterID = ones(numel(t),1); %initialize as noise (ID = 1 being the noise cluster)
pntType = -1*ones(numel(t),1); %initialize as noise

%starting from the long lasting immobilization events as these will be picked up with higher fidelity
%(at smaller searchTime, the chance for false identification is increased due to the less stringent requirements)
for idxT = numel(searchTime):-1:1
    numCluster = max(clusterID) - 1;

    %search for pot. links (observations that potentially stem from the same emitter)
    pntNN = DETRIM_pot_link(x(take),y(take),searchRad);

    %apply forward and backward clustering with DETRIM
    [clusterID(take),pntType(take)] = ...
        DETRIM_fwd_rev_cluster(...
        x(take),y(take),t(take),...
        pntNN,...
        searchTime(idxT),... % current time window size [frames] being tested in the hierarchical search
        searchTime(1),... % absolute minimum duration [frames] for a transient event
        'xFactor',xFactor,...
        'idOffset',numCluster);

    %% start next iteration on the set of unclassified points
    take = find(pntType == -1);
end %for

if verbose
    DETRIM_verbose(x,y,t,clusterID)
end %if
end %fun