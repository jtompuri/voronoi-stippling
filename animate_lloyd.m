function out_file = animate_lloyd(image_path, varargin)
% ANIMATE_LLOYD Create an animation of Lloyd iterations (GIF or MP4).
%
% out_file = animate_lloyd('example.png', 'n_stipples', 8000, 'format','gif')
%
% Name-Value options (subset of stippling.m + animation controls):
%   Core/Init:
%     'n_stipples'              numeric
%     'n_iter'                  numeric (max iterations to animate)
%     'eps'                     numeric (early stop threshold in px)
%     'radius'                  numeric (px, used only for optional circle rendering)
%     'poisson_max_radius_factor' numeric >=1; cap for variable-radius (advanced)
%
%   Density/Sampling:
%     'white_cut'               [0,1] (1 disables suppression)
%     'supersample'             integer >=1 (internal work scale)
%     'supersample_method'      'bicubic'|'bilinear'|'nearest'
%
%   Visualization/Output:
%     'format'                  'gif'|'mp4' (default 'gif')
%     'output_file'             path; default: stipplings/animations/<name>.gif|.mp4
%     'frame_delay'             seconds per frame (gif) (default 0.08 ~ 12.5 fps)
%     'frame_rate'              fps (mp4) (default 20)
%     'show_image'              logical; overlay points on the image (default true)
%     'draw_circles'            logical; draw filled circles (slower) instead of markers (default false)
%     'marker_size'             numeric; scatter marker size if draw_circles=false (default 8)
%     'stipple_color'           color; default StipplerConfig.DEFAULT_STIPPLE_COLOR
%     'background_color'        color; default StipplerConfig.DEFAULT_BACKGROUND_COLOR
%     'verbose'                 logical; default true

    % Parse inputs
    p = inputParser;
    addRequired(p, 'image_path', @(s) ischar(s) || isstring(s));
    addParameter(p, 'n_stipples', StipplerConfig.DEFAULT_N_STIPPLES, StipplerConfig.getValidationFunction('n_stipples'));
    addParameter(p, 'n_iter', StipplerConfig.DEFAULT_N_ITER, StipplerConfig.getValidationFunction('n_iter'));
    addParameter(p, 'eps', StipplerConfig.DEFAULT_EPS, StipplerConfig.getValidationFunction('eps'));
    addParameter(p, 'radius', StipplerConfig.DEFAULT_RADIUS, StipplerConfig.getValidationFunction('radius'));
    addParameter(p, 'poisson_max_radius_factor', StipplerConfig.DEFAULT_POISSON_MAX_RADIUS_FACTOR, StipplerConfig.getValidationFunction('poisson_max_radius_factor'));

    addParameter(p, 'white_cut', StipplerConfig.DEFAULT_WHITE_CUT, StipplerConfig.getValidationFunction('white_cut'));
    addParameter(p, 'supersample', StipplerConfig.DEFAULT_SUPERSAMPLE, StipplerConfig.getValidationFunction('supersample'));
    addParameter(p, 'supersample_method', StipplerConfig.DEFAULT_SUPERSAMPLE_METHOD, StipplerConfig.getValidationFunction('supersample_method'));
    addParameter(p, 'auto_contrast', StipplerConfig.DEFAULT_AUTO_CONTRAST, StipplerConfig.getValidationFunction('auto_contrast'));

    addParameter(p, 'format', 'gif', @(s) any(strcmpi(s, {'gif','mp4'})));
    addParameter(p, 'output_file', '', @(s) ischar(s) || isstring(s));
    addParameter(p, 'frame_delay', 0.08, @(x) isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'frame_rate', 20, @(x) isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'show_image', true, @islogical);
    addParameter(p, 'draw_circles', false, @islogical);
    addParameter(p, 'marker_size', 8, @(x) isnumeric(x) && isscalar(x) && x>0);
    addParameter(p, 'stipple_color', StipplerConfig.DEFAULT_STIPPLE_COLOR, @(x) ischar(x) || (isnumeric(x) && numel(x)==3));
    addParameter(p, 'background_color', StipplerConfig.DEFAULT_BACKGROUND_COLOR, @(x) ischar(x) || (isnumeric(x) && numel(x)==3));
    addParameter(p, 'verbose', true, @islogical);
    parse(p, image_path, varargin{:});

    % Unpack
    image_path = char(p.Results.image_path);
    n_stipples = p.Results.n_stipples;
    n_iter = p.Results.n_iter;
    eps = p.Results.eps;
    radius = p.Results.radius;
    poisson_max_radius_factor = p.Results.poisson_max_radius_factor;
    white_cut = p.Results.white_cut;
    supersample = max(1, round(p.Results.supersample));
    supersample_method = char(p.Results.supersample_method);
    auto_contrast = p.Results.auto_contrast;
    fmt = lower(char(p.Results.format));
    out_file = char(p.Results.output_file);
    delay = p.Results.frame_delay;
    fps = p.Results.frame_rate;
    show_image = p.Results.show_image;
    draw_circles = p.Results.draw_circles;
    msize = p.Results.marker_size;
    stipple_color = p.Results.stipple_color;
    background_color = p.Results.background_color;
    verbose = p.Results.verbose;

    if ~exist(image_path, 'file')
        error('Image file not found: %s', image_path);
    end
    if ~StipplerConfig.isValidImageFormat(image_path)
        [~,~,ext] = fileparts(image_path);
        error('Unsupported image format: %s', ext);
    end

    % Load and preprocess
    img = imread(image_path);
    if size(img,3) == 3, img = rgb2gray(img); end
    img = im2single(img);
    
    % Apply automatic contrast enhancement if enabled
    if auto_contrast
        img = StipplerConfig.enhanceContrast(img, verbose);
    end
    
    [H, W] = size(img);
    rho = 1 - img;
    
    % Use adaptive epsilon if user didn't explicitly provide it
    if any(strcmp('eps', p.UsingDefaults))  % User did NOT provide eps explicitly
        % Use adaptive epsilon based on image size
        eps = StipplerConfig.getAdaptiveEps([H, W]);
        if verbose
            fprintf('[INFO] Using adaptive eps = %.4f for image size %dx%d\n', eps, H, W);
        end
    end
    if supersample > 1
        try
            img_p = imresize(img, supersample, supersample_method);
        catch
            img_p = imresize(img, supersample);
        end
        img_p = min(max(img_p, 0), 1);
        rho_p = 1 - img_p; rho_p = min(max(rho_p, 0), 1);
    else
        img_p = img; rho_p = rho;
    end
    if isfinite(white_cut) && white_cut < 1
        rho_p(img_p >= white_cut) = 0;
    end

    % Output
    if isempty(out_file)
        [~, base, ~] = fileparts(image_path);
        ensure_directory(fullfile('stipplings','animations'));
        out_file = fullfile('stipplings','animations', sprintf('%s_%d_lloyd.%s', base, n_stipples, fmt));
    else
        [out_dir,~,~] = fileparts(out_file);
        if ~isempty(out_dir), ensure_directory(out_dir); end
    end

    % Initialize points
    S = Stippler(verbose);
    % Always variable-radius Poisson-disk initialization
    pts = S.poisson_disk_sample(rho_p, n_stipples, 0, true, poisson_max_radius_factor);

    % Figure setup (headless-safe)
    fig = figure('Visible','off','Color', background_color, 'PaperPositionMode','auto');
    ax = axes('Parent', fig, 'Position', [0 0 1 1]);
    hold(ax, 'on');
    set(ax, 'Color', background_color);
    set(ax, 'YDir','reverse');
    axis(ax, 'equal'); axis(ax, 'off');
    xlim(ax, [0.5, W+0.5]); ylim(ax, [0.5, H+0.5]);

    if show_image
        % Show base image; if supersampled, display original resolution
        imagesc(ax, [1, W], [1, H], img); colormap(ax, gray(256));
    end

    % Artist handles
    scatterH = [];
    patchH = [];

    % Animation writers
    isGif = strcmp(fmt, 'gif');
    if ~isGif
        try
            vw = VideoWriter(out_file, 'MPEG-4'); %#ok<TNMLP>
            vw.FrameRate = fps;
            open(vw);
        catch
            % Fallback if MPEG-4 unsupported
            vw = VideoWriter(out_file, 'Motion JPEG AVI');
            vw.FrameRate = fps;
            open(vw);
        end
    end

    % Iterate and capture
    for it = 1:n_iter
        % Lloyd step on supersampled grid
        labels = S.compute_voronoi_labels(pts, size(rho_p));
        new_pts = S.compute_centroids(labels, rho_p, pts);

        % Measure movement and maybe early stop
        d = hypot(new_pts(:,1)-pts(:,1), new_pts(:,2)-pts(:,2));
        mean_disp = mean(d);
        pts = new_pts;

        % Draw (scaled back if supersampled)
        pts_draw = pts;
        if supersample > 1 && ~isempty(pts_draw), pts_draw = pts_draw ./ supersample; end

        % Clear previous artists efficiently
        if ~isempty(scatterH) && isvalid(scatterH), delete(scatterH); end
        if ~isempty(patchH) && isvalid(patchH), delete(patchH); end

        if draw_circles
            % Filled circle patches (slower for large N)
            segments = max(12, StipplerConfig.DEFAULT_PDF_SEGMENTS);
            theta = linspace(0, 2*pi, segments + 1); theta(end) = [];
            ct = cos(theta); st = sin(theta);
            n = size(pts_draw,1);
            V = zeros(n*segments,2); F = zeros(n,segments);
            for i=1:n
                base = (i-1)*segments + (1:segments);
                V(base,1) = pts_draw(i,1) + radius*ct;
                V(base,2) = pts_draw(i,2) + radius*st;
                F(i,:) = base;
            end
            patchH = patch('Parent', ax, 'Faces', F, 'Vertices', V, ...
                           'FaceColor', stipple_color, 'EdgeColor', 'none');
        else
            scatterH = scatter(ax, pts_draw(:,1), pts_draw(:,2), msize, ...
                               'Marker', '.', 'MarkerEdgeColor', stipple_color); 
        end

        % Iteration text
    ttl = sprintf('Lloyd iteration %d/%d  (mean %s = %.3f px)', it, n_iter, char(916), mean_disp);
        title(ax, ttl, 'Color', [0.2 0.2 0.2], 'FontSize', 10, 'Interpreter','none');

        drawnow;
        frame = getframe(fig);
        [A,map] = rgb2ind(frame.cdata, 256);

        if isGif
            if it == 1
                imwrite(A, map, out_file, 'gif', 'LoopCount', Inf, 'DelayTime', delay);
            else
                imwrite(A, map, out_file, 'gif', 'WriteMode', 'append', 'DelayTime', delay);
            end
        else
            writeVideo(vw, frame);
        end

        if verbose
            fprintf('[ANIM] Iteration %d/%d (mean %s=%.3f)\n', it, n_iter, char(916), mean_disp);
        end
        if mean_disp < eps
            if verbose
                fprintf('[ANIM] Early stop at iteration %d (mean %s < eps)\n', it, char(916));
            end
            break;
        end
    end

    if ~isGif
        close(vw);
    end
    close(fig);
    if verbose
        fprintf('[ANIM] Saved animation: %s\n', out_file);
    end
end
