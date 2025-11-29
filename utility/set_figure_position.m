function coordinates = set_figure_position(ratio, factor, position)
% SET_FIGURE_POSITION Calculates the [left bottom width height] coordinates
% for a figure window based on a desired aspect ratio, a screen size factor,
% and a specified position.
%
%   coordinates = set_figure_position(ratio, factor, position)
%
%   Inputs:
%       ratio:      The desired aspect ratio (Width/Height) of the figure.
%       factor:     A scaling factor (0 < factor <= 1) applied to the screen
%                   dimensions to determine the maximum size of the figure.
%       position:   A string specifying the screen location ('north-west',
%                   'north-east', 'south-east', 'south-west', 'center').
%
%   Output:
%       coordinates: A 4-element vector [left, bottom, width, height] in
%                    pixels, suitable for setting a figure's 'Position' property.
%
%   Copyright (c) 2025 Christian Paolo Richter

% Get the current screen size in pixels: [left bottom width height]
scrSize = get(0, 'ScreenSize');
% Calculate the screen's aspect ratio (Screen Width / Screen Height)
scrRatio = scrSize(3)/scrSize(4);

% --- Calculate Figure Dimensions (Width and Height) ---

% Check if the desired ratio matches the screen's ratio
if  ratio == scrRatio
    % If ratios match, the figure size is simply the screen size multiplied by the factor.
    figWidth = factor*scrSize(3);
    figHeight = factor*scrSize(4);
else
    % If ratios don't match, we assume the figure is initially limited by the screen height.

    % Calculate width based on factor, screen height, and desired ratio
    figWidth = factor*scrSize(3)*ratio/scrRatio;
    % Set figure height based on factor and screen height
    figHeight = factor*scrSize(4);

    % Check if the calculated width exceeds the screen width limit (width limited case)
    if figWidth > scrSize(3)
        % Recalculate dimensions based on the width limit.
        figWidth = factor*scrSize(3); % Width is set by factor*screen width
        % Recalculate height to maintain the desired ratio: H = W / ratio
        figHeight = factor*scrSize(4)/ratio*scrRatio;
    end %if
end %if

% --- Calculate Figure Coordinates (Left and Bottom) ---

% Determine the [left bottom width height] coordinates based on the requested position
switch position
    case 'north-west'
        % Top-left corner: x=1, y=ScreenHeight - FigureHeight
        coordinates = [1 scrSize(4)-figHeight figWidth figHeight];
    case 'north-east'
        % Top-right corner: x=ScreenWidth - FigureWidth, y=ScreenHeight - FigureHeight
        coordinates = [scrSize(3)-figWidth scrSize(4)-figHeight figWidth figHeight];
    case 'south-east'
        % Bottom-right corner: x=ScreenWidth - FigureWidth, y=1
        coordinates = [scrSize(3)-figWidth 1 figWidth figHeight];
    case 'south-west'
        % Bottom-left corner: x=1, y=1
        coordinates = [1 1 figWidth figHeight];
    case 'center'
        % Center: Calculate offset for both x and y
        coordinates = [0.5*(scrSize(3)-figWidth) 0.5*(scrSize(4)-figHeight) figWidth figHeight];
    otherwise
        % Default to center if an unknown position is provided (Good practice)
        coordinates = [0.5*(scrSize(3)-figWidth) 0.5*(scrSize(4)-figHeight) figWidth figHeight];
end %switch

end %fun