function [clusterID,pntType] = ...
    DETRIM_fwd_rev_cluster(x,y,t,pntNN,winT,minT,varargin)
% DETRIM_fwd_rev_cluster performs pure forward and reverse-looking clustering using
% the DBSCAN principle, followed by fusion via maximum core point overlap.
%
% This dual clustering approach is essential to mitigate temporal edge effects
% and ensures the complete reconstruction of a single, fragmented immobile event.
%
%   INPUTS:
%   t: vector; timepoint (image frame) at which the molecule has been localized
%   x: vector; x-position of the molecule [px]
%   y: vector; y-position of the molecule [px]
%   pntNN: cell array of vectors; contains for each point the list indices of all potential links
%   winT: scalar; The current time window size [frames] being tested in the hierarchical search.
%   minT: scalar; The absolute minimum duration [frames] for a transient event.
%
%   OPTIONAL PARAMETERS:
%   idOffset: scalar; (DEFAULT: 0). The starting offset for newly assigned cluster IDs (excluding noise ID 1).
%             Used to prevent cluster ID overlap with previous iterations of the hierarchical search.
%   xFactor: scalar; (DEFAULT: 1.0). The Point Density Correction Factor. Must be between 0 and 1.
%            Used to account for fluorophore blinking/missing data by defining the
%            minimum required number of points for core classification as: (winT * xFactor).
% 
%   OUTPUTS:
%   clusterID: vector; Unique ID for fused clusters (ID=1 is noise).
%   pntType: vector; Classification for fused clusters: 1=Core, 0=Border, -1=Noise.
% 
%   Copyright (c) 2025 Christian Paolo Richter
%   University of Osnabrueck
%
%   modified 18.05.2015: fwd & rev cluster fusion via maximum overlap
%   modified 01.06.2015: make sure only one borderpoint per frame is allocated

%%
ip = inputParser;
ip.KeepUnmatched = true;
addRequired(ip,'x')
addRequired(ip,'y')
addRequired(ip,'t')
addRequired(ip,'pntNN')
addRequired(ip,'winT',@(x)isscalar(x) && x > 0)
addRequired(ip,'minT',@(x)isscalar(x) && x > 0)
addParamValue(ip,'idOffset',0,@(x)isscalar(x)) %#ok
addParamValue(ip,'xFactor',1,@(x)isscalar(x) && x > 0) %#ok
parse(ip,x,y,t,pntNN,winT,minT,varargin{:});

idOffset = ip.Results.idOffset;
xFactor = ip.Results.xFactor;

% Calculate the minimum required number of subsequent localizations (density) 
% required for a point to be a core point.
critScore = xFactor*winT;

%% forward
forwardPntNN = DETRIM_hard_thresh(x,y,t,pntNN,1:winT);
forwardPntScore = cellfun(@numel,forwardPntNN);
%beginn propagation with the densest point in the first frame
[~,idxSort] = sortrows([t max(forwardPntScore)-forwardPntScore]); 
[forwardClusterID,forwardPntType] = ...
    DETRIM_group(t,forwardPntNN,forwardPntScore,critScore,idxSort);

%% reverse
revPntNN = DETRIM_hard_thresh(x,y,t,pntNN,-winT:-1);
revPntScore = cellfun(@numel,revPntNN);
%beginn propagation with the densest point in the last frame
[~,idxSort] = sortrows([max(t)-t max(revPntScore)-revPntScore]); 
[revClusterID,revPntType] = ...
    DETRIM_group(t,revPntNN,revPntScore,critScore,idxSort);

%% fusion of forward & reverse
take = ((forwardPntType == 0) & (revPntType == 1)) | ... %reverse core point
    ((forwardPntType == 1) & (revPntType == 0)) | ... %forward core point
    ((forwardPntType == 1) & (revPntType == 1)); %reverse & forward core point

if not(any(take))
    clusterPair = [];
else
    clusterAsso = [colvec(forwardClusterID(take)) colvec(revClusterID(take))];
    % Use accumarray to count the number of points (overlap) for each (FwdID, RevID) pair
    clusterAsso = accumarray(clusterAsso,1);
    % Find pairs where the overlap is maximal for both the row (Fwd) and column (Rev)
    [clusterPair(:,1),clusterPair(:,2)] = find(...
        bsxfun(@eq,clusterAsso,max(clusterAsso,[],1)) & ...
        transpose(bsxfun(@eq,transpose(clusterAsso),...
        max(transpose(clusterAsso),[],1))) & ...
        clusterAsso > 0);

    [idxDuplicate,I] = find_duplicate(clusterPair(:,1));
    if not(isempty(idxDuplicate))
        for i = 1:numel(idxDuplicate)
            I(find(I(:,i),1,'first'),i) = 0; %retain the first pair
        end %for
        clusterPair(any(I,2),:) = []; %remove all other pairs
    end %if
    [idxDuplicate,I] = find_duplicate(clusterPair(:,2));
    if not(isempty(idxDuplicate))
        for i = 1:numel(idxDuplicate)
            I(find(I(:,i),1,'first'),i) = 0; %retain the first pair
        end %for
        clusterPair(any(I,2),:) = []; %remove all other pairs
    end %if
end

numCluster = size(clusterPair,1);

pntType = -1*ones(size(forwardClusterID));
clusterID = ones(size(forwardClusterID));

cnt = 0;
for idxPair = 1:numCluster
    isCoreFwd = forwardClusterID == clusterPair(idxPair,1) & (forwardPntType == 1);
    isCoreBwd = revClusterID == clusterPair(idxPair,2) & (revPntType == 1);

    %fuse the core points of forward & reverse clustering
    isCore = find(isCoreFwd | isCoreBwd);

    if range(t(isCore)) >= minT
        cnt = cnt + 1;

        numObsT = accumarray(t(isCore),1);
        for idxSimultanObsT = rowvec(find(numObsT > 1))
            isSimultanObsT = find(t(isCore) == idxSimultanObsT);

            pntDist = (x(isCore(isSimultanObsT))-mean(x(isCore))).^2 + ...
                (y(isCore(isSimultanObsT))-mean(y(isCore))).^2;
            isCore(isSimultanObsT(pntDist ~= min(pntDist))) = [];
        end %for

        pntType(isCore) = 1;
        clusterID(isCore) = cnt + idOffset + 1;

        %% identify all borderpoints between the corepoint caps
        isBorderFwd = forwardClusterID == clusterPair(idxPair,1) & (forwardPntType == 0);
        isBorderBwd = revClusterID == clusterPair(idxPair,2) & (revPntType == 0);

        isBorder = find(isBorderFwd | isBorderBwd);
        isBorder = isBorder(t(isBorder) > min(t(isCoreFwd)) & ...
            t(isBorder) < max(t(isCoreBwd)) & ...
            not(ismember(t(isBorder),t(isCore))));

        if not(isempty(isBorder))
            %remove multiple borderpoints within one frame
            numObsT = accumarray(t(isBorder),1);
            %discard all multiple observation (only the NN remains)
            for idxSimultanObsT = rowvec(find(numObsT > 1))
                isSimultanObsT = find(t(isBorder) == idxSimultanObsT);

                %take the border point closest to the cluster center
                pntDist = (x(isBorder(isSimultanObsT))-mean([x(isCore);x(isBorder)])).^2 + ...
                    (y(isBorder(isSimultanObsT))-mean([y(isCore);y(isBorder)])).^2;
                isBorder(isSimultanObsT(pntDist ~= min(pntDist))) = [];
            end %for
        end %if

        pntType(isBorder) = 0;
        clusterID(isBorder) = cnt + idOffset + 1;

    end %if
end %for
end %fun