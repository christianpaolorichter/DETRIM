function pntNN = DETRIM_hard_thresh(x,y,t,pntNN,dt)
% DETRIM_hard_thresh applies the hard temporal filtering required by DETRIM
% to the pre-calculated list of spatial neighbors.
%
% This function reduces the list of potential links (`pntNN`) derived solely
% from spatial proximity (`searchRad`) to those that also reside strictly
% within the expected time window (`dt`) around the query point. It also
% handles the crucial step of selecting the closest point when multiple
% localizations occur within the same time frame and spatial radius, ensuring
% a single, unambiguous link per frame.
%
%   INPUTS:
%   t: vector; timepoint (image frame) at which the molecule has been localized
%   x: vector; x-position of the molecule [px]
%   y: vector; y-position of the molecule [px]
%   pntNN: cell array of vectors; contains for each point the list indices of all potential links
%   dt: vector; defines the time window around each point
%
%   OUTPUTS:
%   pntNN: cell array of vectors; The updated list of neighbors, filtered to only include those that fall within the specified time window.
% 
%   Copyright (c) 2025 Christian Paolo Richter
%   University of Osnabrueck
%
%   modified 19.05.2015: fixed a bug potentially leading to a false time-coordinate for the query point

ip = inputParser;
addRequired(ip,'x')
addRequired(ip,'y')
addRequired(ip,'t')
addRequired(ip,'pntNN',@iscell)
addRequired(ip,'dt')

parse(ip,x,y,t,pntNN,dt);

%%
for pntIdx = numel(pntNN):-1:1
    %select all points that lie within the respective timewindow of the
    %query point
    take = ismember(t(pntNN{pntIdx}),t(pntIdx)+dt);

    pntNN{pntIdx} = pntNN{pntIdx}(take);

    if numel(pntNN{pntIdx}) > 1
        %check if there are multiple observations within 1 frame
        numObsT = accumarray(t(pntNN{pntIdx}),1);

        %discard all multiple observation (only the NN remains)
        for idxSimultanObsT = rowvec(find(numObsT > 1))
            isSimultanObsT = find(t(pntNN{pntIdx}) == idxSimultanObsT);

            pntDist = (x(pntNN{pntIdx}(isSimultanObsT))-x(pntIdx)).^2 + ...
                (y(pntNN{pntIdx}(isSimultanObsT))-y(pntIdx)).^2;

            pntNN{pntIdx}(isSimultanObsT(pntDist ~= min(pntDist))) = [];
        end %for
    end %if
end %for
end %fun