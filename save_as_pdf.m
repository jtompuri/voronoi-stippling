% filepath: save_as_pdf.m
function save_as_pdf(points, sizeWH, filename, radius, segments, stipple_color, background_color)
% SAVE_AS_PDF Save stipples as a vector PDF using circle polygons.
% Draws all dots as filled circular polygons in data units (image pixels),
% then exports as true vector graphics (no rasterization).
%
% Parameters (positional):
%   points            Nx2 [x,y] in image pixel coordinates
%   sizeWH            [width, height] in pixels
%   filename          output .pdf path
%   radius            circle radius in pixels
%   segments          number of polygon segments per circle (default: StipplerConfig.DEFAULT_PDF_SEGMENTS)
%   stipple_color     fill color for dots (default: StipplerConfig.DEFAULT_STIPPLE_COLOR)
%   background_color  page/axes background color (default: StipplerConfig.DEFAULT_BACKGROUND_COLOR)
%
% Notes:
% - Coordinates use top-left origin conventions when combined with YDir='reverse'.
% - A minimal padding (by axis limits) avoids clipping strokes at image edges.
% - Uses exportgraphics with ContentType='vector' for true vector output.

    % Input validation
    if ~isnumeric(points) || size(points, 2) ~= 2
        error('points must be an Nx2 numeric array');
    end
    if ~isnumeric(sizeWH) || length(sizeWH) ~= 2 || any(sizeWH <= 0)
        error('sizeWH must be a 2-element vector with positive values');
    end
    if ~ischar(filename) && ~isstring(filename)
        error('filename must be a string or character array');
    end
    if ~isnumeric(radius) || radius <= 0
        error('radius must be a positive number');
    end

    if nargin < 5 || isempty(segments)
        segments = StipplerConfig.DEFAULT_PDF_SEGMENTS;
    end
    if nargin < 6 || isempty(stipple_color)
        stipple_color = StipplerConfig.DEFAULT_STIPPLE_COLOR;
    end
    if nargin < 7 || isempty(background_color)
        background_color = StipplerConfig.DEFAULT_BACKGROUND_COLOR;
    end
    % No configurable padding; use minimal safe defaults
    
    if ~isnumeric(segments) || segments < StipplerConfig.MIN_PDF_SEGMENTS
        error('segments must be a number >= %d', StipplerConfig.MIN_PDF_SEGMENTS);
    end

    % Ensure output directory exists
    [output_dir, ~, ~] = fileparts(filename);
    ensure_directory(output_dir);

    W = sizeWH(1); H = sizeWH(2);
    n = size(points, 1);
    if n == 0
        warning('No points to render; creating an empty PDF.');
    end

    % Build vertices for all circles in one patch to keep file compact
    theta = linspace(0, 2*pi, segments + 1); 
    theta(end) = []; % drop duplicate
    ct = cos(theta); st = sin(theta);

    % Centers
    cx = points(:,1); cy = points(:,2);
    % For each circle, vertices are (cx + r*ct, cy + r*st)
    % Preallocate vertex array (n*segments x 2)
    V = zeros(n * segments, 2);
    % Faces: each row indexes the segments vertices of one circle
    F = zeros(n, segments);

    for i = 1:n
        base = (i-1)*segments;
        V(base + (1:segments), 1) = cx(i) + radius * ct;
        V(base + (1:segments), 2) = cy(i) + radius * st;
        F(i, :) = base + (1:segments);
    end

    % Create figure sized exactly to image dimensions in pixels
    fig = figure('Visible','off', 'Color', background_color, 'PaperPositionMode','auto');
    
    % Set figure size to exact image dimensions
    try
        set(fig, 'Units', 'pixels');
        set(fig, 'Position', [100, 100, W, H]);
    catch
        % Fallback if pixel sizing fails
    end
    
    % Create axes that fill the entire figure with no margins
    ax = axes(fig, 'Units', 'pixels', 'Position', [0, 0, W, H]);
    hold(ax, 'on');
    set(ax, 'Color', background_color);
    
    % Use exact image canvas (no padding) to preserve original size and aspect ratio
    xlim(ax, [0.5, W+0.5]); ylim(ax, [0.5, H+0.5]);
    set(ax, 'YDir','reverse');         % image-style coordinates
    set(ax, 'DataAspectRatio', [1 1 1]); % 1:1 pixel aspect ratio
    axis(ax, 'off');
    
    % Eliminate all margins and insets (only set writable properties)
    set(ax, 'LooseInset', [0 0 0 0]);
    set(ax, 'Position', [0, 0, W, H]);
    
    % Keep clipping on to maintain exact canvas bounds
    if isprop(ax, 'Clipping'), set(ax, 'Clipping','on'); end
    
    % Draw background rectangle to ensure full canvas is exported (preserves whitespace)
    try
        xl = xlim(ax); yl = ylim(ax);
        bg_rect = rectangle('Parent', ax, 'Position', [xl(1), yl(1), diff(xl), diff(yl)], ...
                           'FaceColor', background_color, 'EdgeColor', 'none');
        % Send background to back
        try, uistack(bg_rect, 'bottom'); catch, end
    catch
        % Best effort; export may crop if this fails
    end

    % Draw as a single black filled patch without edges
    if ~isempty(V)
    patch('Faces', F, 'Vertices', V, 'FaceColor', stipple_color, 'EdgeColor','none', 'Parent', ax);
    % Let axes clipping policy govern overall behavior
    end

    % Export as vector PDF - export the figure to get exact sizing
    try
        if exist('exportgraphics','file')
            exportgraphics(fig, filename, 'ContentType','vector', 'BackgroundColor', background_color);
        else
            print(fig, filename, '-dpdf', '-vector');
        end
    catch ME
        close(fig);
        rethrow(ME);
    end
    close(fig);
    fprintf('[INFO] PDF (vector) saved: %s\n', filename);
end
