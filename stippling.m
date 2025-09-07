function [pdf_file, tsp_file] = stippling(image_path, varargin)
% Weighted Voronoi stippling pipeline
% Usage: [pdf_file, tsp_file] = stippling('image.png', 'n_stipples', 5000, 'radius', 2.0)
% Key options: n_stipples, radius, eps (0=adaptive), supersample

    p = inputParser;
    addRequired(p, 'image_path', @ischar);
    addParameter(p, 'output_basename', '', @ischar);
    addParameter(p, 'n_stipples', StipplerConfig.DEFAULT_N_STIPPLES, StipplerConfig.getValidationFunction('n_stipples'));
    addParameter(p, 'n_iter', StipplerConfig.DEFAULT_N_ITER, StipplerConfig.getValidationFunction('n_iter'));
    addParameter(p, 'radius', StipplerConfig.DEFAULT_RADIUS, StipplerConfig.getValidationFunction('radius'));
    addParameter(p, 'verbose', StipplerConfig.DEFAULT_VERBOSE, StipplerConfig.getValidationFunction('verbose'));
    addParameter(p, 'eps', StipplerConfig.DEFAULT_EPS, StipplerConfig.getValidationFunction('eps'));
    addParameter(p, 'pdf_segments', StipplerConfig.DEFAULT_PDF_SEGMENTS, StipplerConfig.getValidationFunction('pdf_segments'));
    addParameter(p, 'write_pdf', StipplerConfig.DEFAULT_WRITE_PDF, StipplerConfig.getValidationFunction('write_pdf'));
    addParameter(p, 'white_cut', StipplerConfig.DEFAULT_WHITE_CUT, StipplerConfig.getValidationFunction('white_cut'));
    addParameter(p, 'stipple_color', StipplerConfig.DEFAULT_STIPPLE_COLOR, StipplerConfig.getValidationFunction('stipple_color'));
    addParameter(p, 'background_color', StipplerConfig.DEFAULT_BACKGROUND_COLOR, StipplerConfig.getValidationFunction('background_color'));
    addParameter(p, 'supersample', StipplerConfig.DEFAULT_SUPERSAMPLE, StipplerConfig.getValidationFunction('supersample'));
    addParameter(p, 'supersample_method', StipplerConfig.DEFAULT_SUPERSAMPLE_METHOD, StipplerConfig.getValidationFunction('supersample_method'));
    addParameter(p, 'auto_contrast', StipplerConfig.DEFAULT_AUTO_CONTRAST, StipplerConfig.getValidationFunction('auto_contrast'));
    addParameter(p, 'poisson_max_radius_factor', StipplerConfig.DEFAULT_POISSON_MAX_RADIUS_FACTOR, StipplerConfig.getValidationFunction('poisson_max_radius_factor'));
    parse(p, image_path, varargin{:});

    output_basename = p.Results.output_basename;
    n_stipples = p.Results.n_stipples;
    n_iter = p.Results.n_iter;
    radius = p.Results.radius;
    verbose = p.Results.verbose;
    eps = p.Results.eps;
    pdf_segments = p.Results.pdf_segments;
    write_pdf = p.Results.write_pdf;
    stipple_color = p.Results.stipple_color;
    background_color = p.Results.background_color;
    white_cut = p.Results.white_cut;
    supersample = max(1, round(p.Results.supersample));
    supersample_method = char(p.Results.supersample_method);
    auto_contrast = p.Results.auto_contrast;
    poisson_max_radius_factor = p.Results.poisson_max_radius_factor;

    % Input validation
    if ~exist(image_path, 'file')
        error('Image file does not exist: %s', image_path);
    end
    
    % Check if file is a valid image format
    if ~StipplerConfig.isValidImageFormat(image_path)
        [~, ~, ext] = fileparts(image_path);
        error('Unsupported image format: %s. Supported formats: %s', ...
              ext, StipplerConfig.getSupportedFormatsString());
    end
    
    fprintf('[START] Stippling: %s\n', image_path);
    start_time = tic;

    try
        % Ensure output directories
        config = StipplerConfig.getDefault();
        for di = 1:numel(config.output_dirs)
            ensure_directory(config.output_dirs{di});
        end

        % Load & preprocess
        try
            img = imread(image_path);
        catch ME
            error('Failed to read image file "%s": %s', image_path, ME.message);
        end
        if size(img,3) == 3, img = rgb2gray(img); end
        img = im2single(img);
        
        % Apply automatic contrast enhancement if enabled
        if auto_contrast
            img = StipplerConfig.enhanceContrast(img, verbose);
        end
        
    [H, W] = size(img);
    rho = 1 - img;
        
        % Use adaptive epsilon if user didn't explicitly provide it
        if any(strcmp('eps', p.UsingDefaults)) || eps == 0
            eps = StipplerConfig.getAdaptiveEps([H, W]);
            if verbose
                fprintf('[INFO] Using adaptive eps = %.4f for image size %dx%d\n', eps, H, W);
            end
        end
        % Supersample image/density if requested
        if supersample > 1
            try
                img_p = imresize(img, supersample, supersample_method);
            catch
                img_p = imresize(img, supersample); % fallback method default
            end
            % Clamp to [0,1] to avoid resize overshoot
            img_p = min(max(img_p, 0), 1);
            rho_p = 1 - img_p;
            rho_p = min(max(rho_p, 0), 1);
        else
            img_p = img; rho_p = rho;
        end
        % Suppress near-white regions in the density to avoid sampling/relaxation there
        if ~isempty(white_cut) && isfinite(white_cut) && white_cut < 1
            white_mask = img_p >= white_cut;
            rho_p(white_mask) = 0;
            rho_p = min(max(rho_p, 0), 1);
        end
        sizeWH = [W, H];

        if isempty(output_basename)
            [~, input_name, ~] = fileparts(image_path);
            output_basename = sprintf('%s_%d', input_name, n_stipples);
        end

        % Stippling
    stippler = Stippler(verbose);
    % Always variable-radius Poisson-disk initialization
    points = stippler.poisson_disk_sample(rho_p, n_stipples, 0, true, poisson_max_radius_factor);
    [points, iters_done] = stippler.lloyd_relaxation(points, rho_p, n_iter, verbose, eps);
    
    % Scale points back down if supersampled
        if supersample > 1 && ~isempty(points)
            points = points ./ supersample;
        end

        % Save TSP
        tsp_file = fullfile('stipplings','tsp',[output_basename '.tsp']);
        save_as_tsp(points, tsp_file, upper(output_basename));

        % Save PDF
        if write_pdf
            pdf_file = fullfile('stipplings','pdf',[output_basename '.pdf']);
            save_as_pdf(points, sizeWH, pdf_file, radius, pdf_segments, stipple_color, background_color);
        else
            pdf_file = '';
        end

        total_time = toc(start_time);
        fprintf('[RESULT] Iterations: %d\n', iters_done);
        fprintf('[RESULT] Total stipples: %d\n', size(points,1));
        fprintf('[RESULT] Time: %.2f s\n', total_time);
        fprintf('[DONE] Stippling complete.\n');
    catch ME
        fprintf('[ERROR] %s\n', ME.message);
        rethrow(ME);
    end
end
