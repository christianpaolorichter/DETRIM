function [pntNN,pntDist] = DETRIM_pot_link(x,y,searchRad)
% DETRIM_pot_link utilizes a k-d tree search to efficiently establish all
% purely spatial density-reachability links for the DBSCAN-like clustering.
%
% This function determines the initial set of "potential links" or neighbors
% for every localization based solely on the maximum physical distance
% (`searchRad`). This highly optimized pre-filtering step is agnostic to time
% and is essential for speed when dealing with millions of localizations.
%
%   INPUTS:
%   x: vector; x-position of the molecule [px]
%   y: vector; y-position of the molecule [px]
%   searchRad: scalar; defines the max. distance [px] between observations to be considered as potential links
% 
%   OUTPUTS:
%   pntNN: cell array of vectors; For each point, a list of indices (IDs) of all neighbors within searchRad.
%   pntDist: cell array of vectors; For each point, the distances to the corresponding neighbors in pntNN.
%
%   Copyright (c) 2025 Christian Paolo Richter
%   University of Osnabrueck
%
%   modified 08.05.2015: added conversion to uint32 to save 50%

%%
ip = inputParser;
addRequired(ip,'x')
addRequired(ip,'y')
addRequired(ip,'searchRad',@(x)isscalar(x) && x > 0)
parse(ip,x,y,searchRad);

%% calculate the nearest-neighbor relationship
objKdTree = KDTreeSearcher([x,y]);
if nargout == 1
    pntNN = rangesearch(objKdTree,[x,y],searchRad);
else
    [pntNN,pntDist] = rangesearch(objKdTree,[x,y],searchRad);
end %if

if numel(x) <= intmax("uint32")
    pntNN = cellfun(@uint32,pntNN,'un',0); %saves 50% RAM compared to double precision
end %if
end %fun
