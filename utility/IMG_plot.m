function [hImg,hAx,hFig,OUT] = IMG_plot(img,varargin)
% IMG_PLOT Creates a MATLAB figure and axes, then plots a 2D image array
%   using a 'surface' object. It handles coordinate data (xdata/ydata), color
%   scaling (Clim), colormaps, transparency, and a color bar, often used for
%   displaying scientific or processed data with custom settings.
%
%   Inputs:
%       img (matrix): The 2D image data to be plotted (e.g., intensity map).
%       varargin: Optional parameter/value pairs for customizing the plot.
%
%   Outputs:
%       hImg (handle): Handle to the surface object (the image plot).
%       hAx (handle): Handle to the axes object.
%       hFig (handle): Handle to the figure object.
%       OUT (struct): Structure containing handles to other created objects (e.g., color bar).
%
%   Copyright (c) 2014-2025 Christian Paolo Richter
%   University of Osnabrueck
%modified 14.10.2014
%modified 24.11.2014: ip
%modified 13.10.2015

% Initialize the input parser object for handling optional name/value pair arguments.
ip = inputParser;
% Allow arguments not explicitly defined to be stored in ip.Unmatched.
ip.KeepUnmatched = true;

% Define 'img' as a required argument.
addRequired(ip,'img')

% Define optional parameters with default values:
addParamValue(ip,'xdata',[]) % X-coordinates for the image (defaults to 1:width).
addParamValue(ip,'ydata',[]) % Y-coordinates for the image (defaults to 1:height).
addParamValue(ip,'hImg',[])  % Handle to an existing image object to update.
addParamValue(ip,'hAx',[])   % Handle to an existing axes object to plot into.
% Figure size relative to screen size (scalar between 0 and 1).
addParamValue(ip,'FigSize',0.75, @(x)isscalar(x) && x > 0 && x <= 1)
% Saturation limits for color scaling (e.g., [0.01 0.99] for 1st and 99th percentile).
addParamValue(ip,'SatClim',[0.01 0.99])
addParamValue(ip,'FixClim',[]) % Fixed color limits [min max]. Overrides SatClim if used.
addParamValue(ip,'Transparency',[]) % Transparency data (alpha channel).
addParamValue(ip,'ColorMap','gray') % Colormap name (e.g., 'gray', 'jet').
% Color bar label (string). Check that input is a character array.
addParamValue(ip,'ColorBar',[],@(x)ischar(x))
addParamValue(ip,'ColorBarFontSize',20)
addParamValue(ip,'ColorBarColor','k') % Color of the color bar ticks/label.
addParamValue(ip,'ColorBarPosition','left') % Position of the color bar.
addParamValue(ip,'EdgeLinestyle','none') % Linestyle for surface edges (usually 'none' for images).
addParamValue(ip,'EdgeColor','k') % Color of the edges if visible.
addParamValue(ip,'Visible','on') % Figure visibility ('on' or 'off').

% Execute the parsing of the input arguments.
parse(ip,img,varargin{:});

% Extract results from the input parser structure to local variables.
xdata = ip.Results.xdata;
% If xdata was not provided, default to 1:width of the image.
if isempty(xdata)
    xdata = 1:size(img,2);
end %if
ydata = ip.Results.ydata;
% If ydata was not provided, default to 1:height of the image.
if isempty(ydata)
    ydata = 1:size(img,1);
end %if
hImg = ip.Results.hImg;
hAx = ip.Results.hAx;
figSize = ip.Results.FigSize;
satClim = ip.Results.SatClim;
fixClim = ip.Results.FixClim;
colorMap = ip.Results.ColorMap;
colorBar = ip.Results.ColorBar;
colorBarFontSize = ip.Results.ColorBarFontSize;
colorBarColor = ip.Results.ColorBarColor;
colorBarPos = ip.Results.ColorBarPosition;
visible = ip.Results.Visible;
transparency = ip.Results.Transparency;
% If transparency data is not provided, use an opaque matrix of zeros.
if isempty(transparency)
    transparency = zeros(size(img));
end %if

% Initialize the output structure for handles to extra elements (e.g., colorbar).
OUT = [];

%% Plotting Logic
% Check if a handle to an existing image object (hImg) was NOT provided.
if isempty(hImg)
    % If no existing image handle, check if an existing axes handle (hAx) was NOT provided.
    if isempty(hAx)
        %% initialize figure (Create a new figure and axes)

        % Get image dimensions.
        [imgHeight,imgWidth] = size(img);

        % Calculate the figure position based on image aspect ratio and desired relative size.
        % The 'set_figure_position' function (assumed to be available) handles this.
        figPos = set_figure_position(imgWidth/imgHeight, figSize, 'center');

        % Create a new figure with specified properties.
        hFig = figure(...
            'Units','pixels',...
            'Position',figPos,...
            'Color',[1 1 1],... % White background.
            'Visible',visible);
        % Create new axes that fill the entire figure window.
        hAx = axes(...
            'Parent', hFig,...
            'Color',[1 1 1],... % White axes background.
            'CLimMode','auto',...
            'Units','normalized',...
            'Position', [0 0 1 1],... % Axes span 100% of the figure.
            'XTickLabel','',...
            'YTickLabel','',...
            'NextPlot','add',... % Keep existing plots when adding new ones.
            'Box','on');
    else
        % If hAx was provided but hImg was not, get the figure handle from the existing axes.
        hFig = get(hAx,'Parent');
    end %if

    %% Prepare Coordinates and Plot (Plotting using surface)

    % Determine the spacing (pixel width) in the x-direction.
    if numel(xdata) > 1
        dx = xdata(2)-xdata(1);
    else
        % Default spacing if xdata is a single point or empty (handled above).
        dx = 0.5;
    end %if
    % Determine the spacing (pixel height) in the y-direction.
    if numel(ydata) > 1
        dy = ydata(2)-ydata(1);
    else
        % Default spacing if ydata is a single point or empty (handled above).
        dy = 0.5;
    end %if

    % Create the coordinate grids for the surface plot.
    % The coordinates are expanded by one element (xdata(end)+dx, ydata(end)+dy)
    % to define the corners of the grid cells, which is required by surface()
    % when cdata is a matrix of intensity values (MxN array).
    % 'rowvec' (assumed function) ensures xdata/ydata are row vectors for concatenation.
    [xdata,ydata] = meshgrid(...
        [rowvec(xdata) xdata(end)+dx],...
        [rowvec(ydata) ydata(end)+dy]);

    % Create the image using the 'surface' function (more flexible than 'imagesc').
    % The x/y data are shifted by half a pixel (-dx/2, -dy/2) to center the
    % image data ('cdata') over the coordinate points, matching how 'imagesc' works.
    % zdata is set to 0 to ensure a 2D plot.
    hImg = surface(...
        'xdata',xdata-dx/2,...
        'ydata',ydata-dy/2,...
        'zdata',ydata*0,...
        'cdata',double(img)); % Ensure cdata is double precision.

    % Set properties for the newly created surface object (hImg).
    set(hImg,...
        'FaceAlpha','flat',... % Use AlphaData for transparency.
        'AlphaData',1-transparency,... % Alpha is 1 - Transparency (1=Opaque, 0=Transparent).
        'alphadatamapping','none',...
        'CDatamapping','scaled',... % Map CData (image values) to color limits.
        'linestyle',ip.Results.EdgeLinestyle,...
        'EdgeColor',ip.Results.EdgeColor,...
        'Parent',hAx); % Link the image to the axes.

else
    % This block executes if a handle to an existing image object (hImg) was provided,
    % meaning the function is UPDATING an existing plot.

    % Retrieve handles for the parent axes and figure.
    hAx = get(hImg,'Parent');
    hFig = get(hAx,'Parent');

    % Recalculate dx and dy spacings (same logic as above).
    if numel(xdata) > 1
        dx = xdata(2)-xdata(1);
    else
        dx = 0.5;
    end %if
    if numel(ydata) > 1
        dy = ydata(2)-ydata(1);
    else
        dy = 0.5;
    end %if

    % Recalculate coordinate grids for the updated surface plot (same logic as above).
    [xdata,ydata] = meshgrid(...
        [rowvec(xdata) xdata(end)+dx],...
        [rowvec(ydata) ydata(end)+dy]);

    % Update properties of the existing surface object (hImg).
    set(hImg,...
        'xdata',xdata-dx/2,...
        'ydata',ydata-dy/2,...
        'zdata',ydata*0,...
        'cdata',img,... % Update image data.
        'FaceAlpha','flat',...
        'AlphaData',1-transparency,... % Update transparency.
        'alphadatamapping','none',...
        'CDatamapping','scaled',...
        'linestyle',ip.Results.EdgeLinestyle,...
        'EdgeColor',ip.Results.EdgeColor);
end %if

% --- Common Post-Plotting Operations (Executed for both new and updated plots) ---

% Set the colormap for the axes object.
colormap(hAx,colorMap)

% Set the color limits (caxis) which define how image values are mapped to colors.
try
    if fixClim
        % Option 1: Use explicitly fixed color limits.
        caxis(hAx,fixClim)
    else
        % Option 2: Use saturated color limits by calculating quantiles of the image data.
        % The 'nanquantile' function (assumed to be available) handles NaNs.
        % This effectively ignores the highest/lowest 'SatClim' fraction of data
        % to prevent outliers from skewing the color mapping.
        caxis(hAx,nanquantile(img,satClim));
    end %if
end %try (Used to suppress errors if quantile or caxis fails)

% Handle Color Bar creation/update.
if not(isempty(colorBar))
    % Check requested position and create the color bar object.
    switch colorBarPos
        case 'left'
            OUT.hCbar = colorbar(hAx,'Location','west');
        case 'bottom'
            OUT.hCbar = colorbar(hAx,'Location','South');
        otherwise
            OUT.hCbar = colorbar(hAx,'Location',colorBarPos);
    end %switch

    % Customize the appearance of the color bar ticks and label.
    set(OUT.hCbar,'xcolor',colorBarColor,'ycolor',colorBarColor,...
        'fontsize',colorBarFontSize,'fontweight','bold')

    % Set the Y-label text and color on the color bar.
    set(get(OUT.hCbar,'YLabel'),'String',colorBar,'Color',colorBarColor)
end %if

% Ensure the figure is rendered immediately (important for scripts or functions).
drawnow expose
end %fun