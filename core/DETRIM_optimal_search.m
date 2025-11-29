function [searchRad,searchTime] = DETRIM_optimal_search(D,locPrec,varargin)
% DETRIM_OPTIMAL_SEARCH determines the minimum observation time window (searchTime)
%   and corresponding optimal search radius (searchRad) required for the DETRIM
%   algorithm to achieve a user-defined False Positive Rate (typically 1%),
%   given the physical properties of the system (Diffusion coefficient and
%   Localization precision).
%
%   INPUTS:
%   D: scalar; Diffusion coefficient [px^2/frame]. (A measure of particle mobility.)
%   locPrec: scalar; Localization precision [px]. (A measure of measurement noise.)
%   varargin: Optional parameter/value pairs, primarily used to set 'targetFPR'.
%
%   OPTIONAL PARAMETERS:
%   targetFPR: scalar; (DEFAULT: 0.01). The maximum acceptable False Positive Rate (FPR)
%              for the DETRIM detection method. FPR is the probability of incorrectly
%              classifying a mobile/noise point as an immobile one.
%
%   OUTPUTS:
%   searchRad: scalar; The optimal spatial search radius [px] determined by the helper
%              function for the final searchTime.
%   searchTime: integer; The minimum number of sequential frames (time window)
%               required to satisfy the target FPR.
%
%   Copyright (c) 2025 Christian Paolo Richter
%   University of Osnabrueck

% Initialize the input parser object for flexible argument handling.
ip = inputParser;
ip.KeepUnmatched = true;
addRequired(ip,'D')
addRequired(ip,'locPrec')
addParamValue(ip,'targetFPR',0.01)
parse(ip,D,locPrec,varargin{:});

% Initialize the temporal search window size.
searchTime = 1;

% Initialize the False Positive Rate (FPR) to infinity to guarantee that the optimization loop begins.
FPR = inf;

% Main optimization loop: Iteratively increases 'searchTime' until the calculated FPR is
% less than or equal to the target FPR.
while FPR > ip.Results.targetFPR

    % Increment the search time window by one frame. Increasing time generally improves detection reliability.
    searchTime = searchTime + 1;

    % simulate trajectories (based on D, locPrec) and
    % evaluate the DETRIM performance for the current 'searchTime'.
    [FPR,searchRad] = sim_and_apply_DETRIM(...
        D,locPrec,searchTime);

end %while
end %fun