function [clusterID,pntType,numCluster,clusterSize] = DETRIM_group(t,pntNN,pntScore,critScore,pntOrder)
% DETRIM_group implements the specialized DBSCAN-like clustering required
% for DETRIM's sequential processing.
%
% This function assigns points to a cluster (or labels them as noise) by
% propagating clusters via density-reachability. Crucially, it ensures
% that clusters only propagate through Core Points (`pntType == 1`) to
% maintain cluster integrity.
%
%   INPUTS:
%   t: vector; timepoint (image frame) at which the molecule has been localized
%   pntNN: cell array of vectors; contains for each point the list indices of all potential links (neighbors)
%   pntScore: vector; density score (number of neighbors in winTime) for each point
%   critScore: scalar; threshold for core point classification (MinPoints required)
%   pntOrder: vector; The order in which points are processed (sorted by time then density).
%
%   OUTPUTS:
%   clusterID: vector; Local cluster identifier (1, 2, 3...) for points classified as Core/Border.
%   pntType: vector; Classification: 1=Core Point, 0=Border Point, -1=Noise/Mobile Point.
%   numCluster: scalar; Total number of clusters found (excluding noise).
%   clusterSize: vector; Size (number of points) of each cluster.
%
%   Copyright (c) 2025 Christian Paolo Richter
%   University of Osnabrueck
%
%   modified 18.05.2015: fixed a bug that allowed the clusters to propagate through border points

ip = inputParser;
addRequired(ip,'t')
addRequired(ip,'pntNN')
addRequired(ip,'pntScore')
addRequired(ip,'critScore',@(x)isscalar(x) && x > 0)
addRequired(ip,'pntOrder')
parse(ip,t,pntNN,pntScore,critScore,pntOrder);

%%
numPnts = numel(pntScore);

pntType = (pntScore >= critScore); %is core point

clusterID = nan(numPnts,1); %unclassified
clusterIdx = 0;
for pntIdx = rowvec(pntOrder)
    if pntType(pntIdx) == 1 && isnan(clusterID(pntIdx))
        %initialize new cluster
        clusterIdx = clusterIdx+1;

        %point is put into respective cluster
        clusterID(pntIdx) = clusterIdx;

        makeConn = pntNN{pntIdx};
        makeConn(not(isnan(clusterID(makeConn)))) = []; %only use those points not already classified

        clusterID(makeConn) = clusterIdx;
        while any(makeConn)
            makeConn(not(pntType(makeConn) == 1)) = []; %make sure clusters propagate only through core points

            %find all connected points (density-reacheability)
            makeConn = setdiff(horzcat(pntNN{makeConn}),makeConn);
            makeConn(not(isnan(clusterID(makeConn)))) = []; %only use those points not already classified
            
            % Block propagation through time points already in the current cluster
            makeConn(ismember(t(makeConn),...
                t(clusterID == clusterIdx))) = [];

            %use only the strongest connection for propagation
            [~,take] = max(bsxfun(@times,...
                double(bsxfun(@eq,unique(t(makeConn)),...
                rowvec(t(makeConn)))),...
                rowvec(pntScore(makeConn))),[],2);
            makeConn = makeConn(take);
            clusterID(makeConn) = clusterIdx;
        end %while
    end %if
end %for
clusterID = clusterID + 1;
isNoise = isnan(clusterID);
clusterID(isNoise) = 1; %assign noise points to the same cluster

pntType = double(pntType);
pntType(isNoise) = -1; %is noise point

numCluster = clusterIdx; %(= #, noise cluster excluded)
clusterSize = accumarray(clusterID,1);
end %fun