function out_file = visualize_tsp_tour(tsp_input, varargin)
% VISUALIZE_TSP_TOUR Render a TSP tour as a closed path (top-left origin)
%
% Minimal, tour-only visualization mirroring visualize.py behavior.
%
% Usage:
%   visualize_tsp_tour('example-1024px_100')
%   visualize_tsp_tour('example-1024px_100', 'line_width', 2, 'join_style','round')
%   visualize_tsp_tour('stipplings/tsp/example-1024px_100.tsp', 'pdf_output','visualizations/pdf/out.pdf')
%
% Inputs:
%   tsp_input  - Basename (without extension) OR path to a .tsp file.
%
% Name-Value Parameters:
%   'tsp_path'       - Full path to .tsp (overrides tsp_input)
%   'tour_path'      - Full path to .tour (if omitted, inferred from TSP path)
%   'output_file'    - Image output (e.g., .png). Optional.
%   'pdf_output'     - Vector PDF output path. Optional.
%   'line_width'     - Tour line width (default: 2.0)
%   'line_color'     - Tour line color (default: 'k')
%   'background_color' - Figure and page background color (default: 'white')
%   'join_style'     - 'miter' (default) | 'round' | 'bevel' | 'chamfer'
%   'miter_limit'    - Miter limit for 'miter' joins (default: 10)
%   'verbose'        - Print progress (default: true)
%   'use_spline'     - Render a curved path through the tour (default: false)
%   'spline_method'  - 'cscvn' (Spline Toolbox) | 'pchip' (fallback) (default: 'cscvn')
%   'spline_samples' - Number of samples along the curve (default: ~2x tour length, capped)

    % Parse inputs
    p = inputParser;
    addRequired(p, 'tsp_input', @(s) ischar(s) || isstring(s));
    addParameter(p, 'tsp_path', '', @(s) ischar(s) || isstring(s));
    addParameter(p, 'tour_path', '', @(s) ischar(s) || isstring(s));
    addParameter(p, 'output_file', '', @(s) ischar(s) || isstring(s));
    addParameter(p, 'pdf_output', '', @(s) ischar(s) || isstring(s));
    addParameter(p, 'line_width', 2.0, @(x) isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'line_color', StipplerConfig.DEFAULT_STIPPLE_COLOR, @(x) ischar(x) || (isnumeric(x) && numel(x)==3));
    addParameter(p, 'background_color', StipplerConfig.DEFAULT_BACKGROUND_COLOR, @(x) ischar(x) || (isnumeric(x) && numel(x)==3));
    addParameter(p, 'join_style', 'miter', @ischar);
    addParameter(p, 'miter_limit', 10, @(x) isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'verbose', true, @islogical);
    addParameter(p, 'use_spline', false, @islogical);
    addParameter(p, 'spline_method', 'cscvn', @ischar);
    addParameter(p, 'spline_samples', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>10));
    parse(p, tsp_input, varargin{:});

    tsp_path = strtrim(char(p.Results.tsp_path));
    tour_path = strtrim(char(p.Results.tour_path));
    output_file = strtrim(char(p.Results.output_file));
    pdf_output = strtrim(char(p.Results.pdf_output));
    line_width = p.Results.line_width;
    line_color = p.Results.line_color;
    background_color = p.Results.background_color;
    join_style = lower(strtrim(p.Results.join_style));
    if strcmp(join_style,'chamfer'), join_style = 'bevel'; end
    miter_limit = p.Results.miter_limit;
    verbose = p.Results.verbose;
    use_spline = p.Results.use_spline;
    spline_method = lower(strtrim(p.Results.spline_method));
    spline_samples = p.Results.spline_samples;

    % Resolve TSP path: explicit arg wins; else from tsp_input
    if isempty(tsp_path)
        tsp_input = char(p.Results.tsp_input);
        if endsWith(tsp_input, '.tsp') || contains(tsp_input, filesep)
            tsp_path = tsp_input;
        else
            % Assume project layout: stipplings/tsp/<basename>.tsp
            tsp_path = fullfile(pwd, 'stipplings','tsp',[tsp_input '.tsp']);
        end
    end
    if ~exist(tsp_path, 'file')
        error('TSP file not found: %s', tsp_path);
    end

    % Resolve tour path: explicit arg or infer next to TSP
    if isempty(tour_path)
        [tp_dir, tp_name, ~] = fileparts(tsp_path);
        cand = fullfile(tp_dir, [tp_name '.tour']);
        if exist(cand, 'file')
            tour_path = cand;
        else
            error('Tour file not found and could not be inferred from TSP path.');
        end
    end
    if ~exist(tour_path, 'file')
        error('Tour file not found: %s', tour_path);
    end

    if verbose
        fprintf('[INFO] Visualizing tour for: %s\n', tsp_path);
    end

    % Load data
    coords = read_tsp_coords(tsp_path);
    [tour0, tour_header_count] = read_linkern_tour(tour_path); % 0-based
    tour = tour0 + 1; % to 1-based
    if numel(tour) >= 2 && tour(1) == tour(end)
        tour = tour(1:end-1);
    end
    % Validate indices are within coordinate range
    ncoords = size(coords,1);
    % Validate that tour header count matches TSP coordinates
    if ~isnan(tour_header_count) && tour_header_count ~= ncoords
        error(['Tour/TSP mismatch: TSP has %d coords but the tour header says %d. ', ...
               'This often happens if you changed the point set (e.g., 99%% white suppression or duplicate removal) ', ...
               'after generating the tour. Regenerate the tour for the filtered TSP.'], ncoords, tour_header_count);
    end
    if isempty(tour) || any(~isfinite(tour)) || any(tour < 1 | tour > ncoords)
        tmin = min(tour); tmax = max(tour);
        error('Tour has invalid indices after 1-based conversion (min=%d, max=%d, ncoords=%d).', tmin, tmax, ncoords);
    end
    if verbose
        fprintf('[INFO] Loaded %d coords, %d tour nodes\n', size(coords,1), numel(tour));
    end

    % Determine canvas size from coordinate bounds
    if isempty(coords)
        W = 100; H = 100; % fallback
    else
        minX = min(coords(:,1)); maxX = max(coords(:,1));
        minY = min(coords(:,2)); maxY = max(coords(:,2));
        W = ceil(maxX - minX + 1);
        H = ceil(maxY - minY + 1);
        % Ensure minimum size
        W = max(W, 100); H = max(H, 100);
    end

    % Prepare figure sized to canvas dimensions
    fig = figure('Visible','off','Color', background_color, 'PaperPositionMode','auto');
    
    % Set figure size to canvas dimensions
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

    % Set coordinate bounds with padding to avoid clipping tour lines
    if isempty(coords)
        xlim(ax, [0.5, W+0.5]); ylim(ax, [0.5, H+0.5]);
    else
        % Add padding based on line width to prevent clipping
        pad = max(line_width * 2, 5); % ensure sufficient padding for thick lines
        xlim(ax, [min(coords(:,1))-pad, max(coords(:,1))+pad]); 
        ylim(ax, [min(coords(:,2))-pad, max(coords(:,2))+pad]);
    end
    set(ax, 'YDir','reverse');         % image-style coordinates
    set(ax, 'DataAspectRatio', [1 1 1]); % 1:1 pixel aspect ratio
    axis(ax, 'off');
    
    % Eliminate all margins and insets
    set(ax, 'LooseInset', [0 0 0 0]);
    set(ax, 'Position', [0, 0, W, H]);
    
    % Disable clipping to avoid cutting off tour segments
    if isprop(ax, 'Clipping'), set(ax, 'Clipping','off'); end
    
    % Draw background rectangle based on coordinate bounds (not padded limits)
    try
        if isempty(coords)
            bg_x = 0.5; bg_y = 0.5; bg_w = W; bg_h = H;
        else
            bg_x = min(coords(:,1)) - 0.5;
            bg_y = min(coords(:,2)) - 0.5;
            bg_w = max(coords(:,1)) - min(coords(:,1)) + 1;
            bg_h = max(coords(:,2)) - min(coords(:,2)) + 1;
        end
        bg_rect = rectangle('Parent', ax, 'Position', [bg_x, bg_y, bg_w, bg_h], ...
                           'FaceColor', background_color, 'EdgeColor', 'none');
        % Send background to back
        try, uistack(bg_rect, 'bottom'); catch, end
    catch
        % Best effort; export may crop if this fails
    end

    % Build path: straight polyline or spline curve
    tour_x = coords(tour,1); tour_y = coords(tour,2);
    % Ensure explicit closure
    if tour_x(1) ~= tour_x(end) || tour_y(1) ~= tour_y(end)
        tour_x(end+1) = tour_x(1);
        tour_y(end+1) = tour_y(1);
    end

    if use_spline
        % Decide samples: proportional to tour length but capped
        if isempty(spline_samples)
            nbase = numel(tour_x);
            spline_samples = min(30000, max(2000, 2*nbase));
        end
        [cx, cy] = compute_curve(tour_x, tour_y, spline_method, spline_samples);
        h = line(ax, 'XData', cx, 'YData', cy, 'Color', line_color, 'LineWidth', line_width); 
    else
        h = patch('Parent', ax, 'XData', tour_x, 'YData', tour_y, ...
                  'FaceColor','none','EdgeColor', line_color, 'LineWidth', line_width, ...
                  'LineJoin', join_style);
        if strcmp(join_style,'miter')
            try set(h, 'MiterLimit', miter_limit); catch, end
        end
    end
    try set(h, 'Clipping','off'); catch, end

    % Canvas is already sized correctly; no need to adjust limits

    % Export: image and/or PDF; default to project PDF if none provided
    out_file = '';
    if isempty(output_file) && isempty(pdf_output)
        % default PDF under stipplings/pdf/<basename>_tour.pdf
        [~, base_name, ~] = fileparts(tsp_path);
        ensure_directory(fullfile(pwd, 'stipplings','pdf'));
        pdf_output = fullfile(pwd, 'stipplings','pdf',[base_name '_tour.pdf']);
    end

    try
        if ~isempty(output_file)
            % Raster image
            [out_dir,~,~] = fileparts(output_file);
            if ~isempty(out_dir), ensure_directory(out_dir); end
            if exist('exportgraphics','file')
                exportgraphics(ax, output_file, 'ContentType','image', 'BackgroundColor', background_color, 'Resolution', 300, 'Padding',0);
            else
                print(fig, output_file, '-dpng', '-r300');
            end
            out_file = output_file;
            if verbose, fprintf('[INFO] Saved image: %s\n', output_file); end
        end
        if ~isempty(pdf_output)
            [out_dir,~,~] = fileparts(pdf_output);
            if ~isempty(out_dir), ensure_directory(out_dir); end
            if exist('exportgraphics','file')
                exportgraphics(fig, pdf_output, 'ContentType','vector', 'BackgroundColor', background_color);
            else
                print(fig, pdf_output, '-dpdf', '-vector');
            end
            out_file = pdf_output;
            if verbose, fprintf('[INFO] Saved PDF: %s\n', pdf_output); end
        end
    catch ME
        close(fig); rethrow(ME);
    end
    close(fig);
end

function coordinates = read_tsp_coords(tsp_file)
% Read coordinates from TSP file
    coordinates = [];
    fid = fopen(tsp_file, 'r');
    if fid == -1
        error('Cannot open TSP file: %s', tsp_file);
    end
    
    try
        % Skip header until NODE_COORD_SECTION (case-insensitive)
        while ~feof(fid)
            line = fgetl(fid);
            if ~ischar(line), break; end
            if contains(lower(strtrim(line)), 'node_coord_section')
                break;
            end
        end

        % Read coordinates; tolerate extra spaces/tabs/comments
        numericLine = @(s) ~isempty(regexp(s, '^\s*\d+\s+[-+]?\d*\.?\d+(e[-+]?\d+)?\s+[-+]?\d*\.?\d+(e[-+]?\d+)?\s*$', 'once'));
        while ~feof(fid)
            line = fgetl(fid);
            if ~ischar(line), break; end
            t = strtrim(line);
            if isempty(t), continue; end
            if contains(lower(t), 'eof'), break; end
            if ~numericLine(t), continue; end
            parts = sscanf(t, '%d %f %f');
            if numel(parts) == 3
                coordinates = [coordinates; parts(2), parts(3)]; %#ok<AGROW>
            end
        end

        fclose(fid);
    catch ME
        fclose(fid);
        rethrow(ME);
    end
end

function [tour, node_count] = read_linkern_tour(tour_file)
% Read a Linkern tour file (first line is count, then pairs: "curr next")
    fid = fopen(tour_file, 'r');
    if fid == -1, error('Cannot open tour file: %s', tour_file); end
    cleaner = onCleanup(@() fclose(fid));
    % First line is node count (Linkern format). Keep it for validation.
    header = fgetl(fid);
    node_count = NaN;
    if ischar(header)
        hc = str2double(strtrim(header));
        if isfinite(hc) && hc > 0
            node_count = hc;
        end
    end
    % Read the first pair
    line = fgetl(fid);
    if ~ischar(line), error('Empty tour file: %s', tour_file); end
    nums = sscanf(strtrim(line), '%d %d');
    if numel(nums) < 2, error('Invalid tour file line: %s', line); end
    current = nums(1);
    tour = current; % 0-based
    % Append next nodes from second column of each line
    while true
        if ~ischar(line), break; end
        nums = sscanf(strtrim(line), '%d %d');
        if numel(nums) >= 2
            nextn = nums(2);
            tour(end+1) = nextn; %#ok<AGROW>
        end
        line = fgetl(fid);
    end
end

function [cx, cy] = compute_curve(x, y, method, n_samples)
% Compute a closed smooth curve through/near points.
% method: 'cscvn' (interpolating cubic spline if available) or 'pchip' fallback.
    x = x(:)'; y = y(:)';
    P = [x; y];
    cx = []; cy = [];
    if strcmp(method, 'cscvn') && exist('cscvn','file')
        try
            fn = cscvn(P);
            t = linspace(fn.breaks(1), fn.breaks(end), n_samples);
            C = fnval(fn, t);
            cx = C(1,:); cy = C(2,:);
        catch
        end
    end
    if isempty(cx)
        % Fallback: parametrize by cumulative chord length and use PCHIP
        s = [0, cumsum(hypot(diff(x), diff(y)))];
        % Ensure strictly increasing for pchip
        [s, iu] = unique(s, 'stable');
        xu = x(iu); yu = y(iu);
        tt = linspace(0, s(end), n_samples);
        cx = pchip(s, xu, tt);
        cy = pchip(s, yu, tt);
    % pchip fallback used
    end
    % Ensure closed-curve exact closure at ends
    if cx(1) ~= cx(end) || cy(1) ~= cy(end)
        cx(end) = cx(1); cy(end) = cy(1);
    end
end
