function DETRIM_verbose(x,y,t,clusterID)
% DETRIM_VERBOSE visualizes the results of the DETRIM cluster result
%
%   INPUTS:
%   t: vector; timepoint (image frame) at which the molecule has been localized
%   x: vector; x-position of the molecule
%   y: vector; y-position of the molecule
%
%   Copyright (c) 2025 Christian Paolo Richter
%   University of Osnabrueck

%% visualization
% Create a new figure window with a white background and size it for three panels
hFig = figure('color','w','Units','pixels',...
    'Position',set_figure_position(3,0.5,'center'),'visible','on');

% --- Panel 1: 2D Projection (x vs y) ---
% Define the first axes for the 2D plot (left third of the figure)
hAx = axes('units','normalized','outerposition',[0 0 1/3 1],'box','on');
hold on

% Determine the unique cluster IDs (ID=1 is noise)
numCluster = max(clusterID)-1; % Count of valid clusters (IDs > 1)

% plot noise points
take = (clusterID == 1);
plot3(x(take),y(take),t(take),'k.','markersize',9)

if numCluster > 0
    % Generate distinct colors for each valid cluster using HSV to RGB conversion.
    % randperm is used to shuffle the colors, ensuring adjacent clusters have different hues.
    clusterColor = hsv2rgb([colvec(randperm(numCluster))/...
        numCluster,ones(numCluster,2)]);
    
    % Loop through all valid cluster IDs (starting from the 2nd unique ID, index 2)
    for idxCluster = 2:max(clusterID)
        take = (clusterID == idxCluster);
        
        % Plot core points for the current cluster (colored circle marker)
        plot3(x(take),y(take),t(take),'.',...
            'color',clusterColor(idxCluster-1,:),'markersize',10)        
    end %for
end %if

% Set axis labels and aspect ratio for the 2D projection
xlabel(hAx(1),'x [px]');
ylabel(hAx(1),'y [px]');
axis image ij % Sets aspect ratio for spatial coordinates

%%
% --- Panel 2: 3D Spatio-Temporal View (x vs y vs t) ---
% Duplicate the first axes object to preserve the plotted data
hAx(2) = copyobj(hAx,hFig); 
% Reposition the copied axis to the middle third of the figure
set(hAx(2),'outerposition',[0.375 0.05 1/3*0.9 0.9],...
    'DataAspectRatio',[range(x)/range(t) range(y)/range(t) 1])
% Change the view to 3D
view(hAx(2),3)
axis(hAx(2),'vis3d')
% Set axis labels for the 3D plot
xlabel(hAx(2),'x [px]');
ylabel(hAx(2),'y [px]');
zlabel(hAx(2),'Time [frame]')

%%
% --- Panel 3: Cluster Duration Histogram ---
% Define the third axes for the histogram (right third of the figure)
hAx(3) = axes('units','normalized','outerposition',[2/3 0 1/3 1],'box','on');

tauObs = []; % Initialize array to store observed cluster durations

% Iterate backwards through cluster IDs (2 is the first valid cluster)
for idxCluster = max(clusterID):-1:2
    take = clusterID == idxCluster;
    % Calculate the observed duration (range in frames) for each cluster
    tauObs(idxCluster,1) = range(t(take));
end %for
tauObs(1) = []; % remove noise cluster duration (index 1 corresponds to ID=1)

% Define histogram bins (duration length)
bins = 1:max(tauObs);
% Calculate frequency (N) of each observed duration (1 frame, 2 frames, etc.)
N = accumarray(tauObs,1); 

% Plot histogram of observed durations using scatter markers
stem(bins,N,'filled','color','k','MarkerFaceColor','k',...
     'MarkerEdgeColor','k','marker','.')

% Set histogram labels and title
xlabel(hAx(3),'Observed Immobilization Duration [frame]');
ylabel(hAx(3),'Frequency (N)');
title(hAx(3),sprintf('Total # detected Immobilization Events = %d',sum(N)))
box on
grid on
axis tight
end