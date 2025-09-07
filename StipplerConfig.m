classdef StipplerConfig
    properties (Constant)
        % Algorithm defaults
        DEFAULT_N_STIPPLES = 5000;
        DEFAULT_N_ITER = 50;
        DEFAULT_RADIUS = 2.0;
        DEFAULT_EPS = 0;                % 0 = adaptive
        DEFAULT_WHITE_CUT = 1.0;
        DEFAULT_PDF_SEGMENTS = 24;
        DEFAULT_WRITE_PDF = true;
        DEFAULT_VERBOSE = true;

        % Processing
        DEFAULT_SUPERSAMPLE = 1;
        DEFAULT_SUPERSAMPLE_METHOD = 'bicubic';
        
        % Auto-contrast
        DEFAULT_AUTO_CONTRAST = false;
        CONTRAST_THRESHOLD = 0.3;
        CONTRAST_ENHANCEMENT_FACTOR = 1.5;
        
        DEFAULT_POISSON_MAX_RADIUS_FACTOR = 4.0;
            
        % Visual
        DEFAULT_STIPPLE_COLOR = 'black';
        DEFAULT_BACKGROUND_COLOR = 'white';
        
        % Limits
        MIN_STIPPLES = 1;
        MAX_STIPPLES = 100000;
        MIN_ITERATIONS = 1;
        MAX_ITERATIONS = 500;
        MIN_RADIUS = 0.1;
        MAX_RADIUS = 50.0;
        MIN_PDF_SEGMENTS = 6;
        MAX_PDF_SEGMENTS = 360;
        MIN_EPS = 1e-6;
        MAX_EPS = 10.0;
        
        % Internal
        MAX_ATTEMPTS_MULTIPLIER = 20;
        FALLBACK_EPS = 0.01;
        BASE_EPS = 0.015;               % 512x512 reference
        EPS_SCALE_FACTOR = 0.0025;
        
        % Files
        SUPPORTED_IMAGE_FORMATS = {'.jpg', '.jpeg', '.png', '.bmp', '.tiff', '.tif', '.gif'};
        OUTPUT_BASE_DIR = 'stipplings';
        OUTPUT_PDF_DIR = 'stipplings/pdf';
        OUTPUT_TSP_DIR = 'stipplings/tsp';
    end
    
    methods (Static)
        function config = getDefault()
            config = struct();
            config.n_stipples = StipplerConfig.DEFAULT_N_STIPPLES;
            config.n_iter = StipplerConfig.DEFAULT_N_ITER;
            config.radius = StipplerConfig.DEFAULT_RADIUS;
            config.eps = StipplerConfig.DEFAULT_EPS;
            config.pdf_segments = StipplerConfig.DEFAULT_PDF_SEGMENTS;
            config.verbose = StipplerConfig.DEFAULT_VERBOSE;
            config.write_pdf = StipplerConfig.DEFAULT_WRITE_PDF;
            config.stipple_color = StipplerConfig.DEFAULT_STIPPLE_COLOR;
            config.background_color = StipplerConfig.DEFAULT_BACKGROUND_COLOR;
            config.white_cut = StipplerConfig.DEFAULT_WHITE_CUT;
            config.supersample = StipplerConfig.DEFAULT_SUPERSAMPLE;
            config.supersample_method = StipplerConfig.DEFAULT_SUPERSAMPLE_METHOD;
            config.auto_contrast = StipplerConfig.DEFAULT_AUTO_CONTRAST;
            config.poisson_max_radius_factor = StipplerConfig.DEFAULT_POISSON_MAX_RADIUS_FACTOR;
            config.output_dirs = {
                StipplerConfig.OUTPUT_BASE_DIR,
                StipplerConfig.OUTPUT_PDF_DIR,
                StipplerConfig.OUTPUT_TSP_DIR
            };
        end
        
        function validator = getValidationFunction(param_name)
            switch lower(param_name)
                case 'n_stipples'
                    validator = @(x) isnumeric(x) && x >= StipplerConfig.MIN_STIPPLES && x <= StipplerConfig.MAX_STIPPLES;
                case 'n_iter'
                    validator = @(x) isnumeric(x) && x >= StipplerConfig.MIN_ITERATIONS && x <= StipplerConfig.MAX_ITERATIONS;
                case 'radius'
                    validator = @(x) isnumeric(x) && x >= StipplerConfig.MIN_RADIUS && x <= StipplerConfig.MAX_RADIUS;
                case 'eps'
                    validator = @(x) isnumeric(x) && x >= 0 && x <= StipplerConfig.MAX_EPS;
                case 'pdf_segments'
                    validator = @(x) isnumeric(x) && x >= StipplerConfig.MIN_PDF_SEGMENTS && x <= StipplerConfig.MAX_PDF_SEGMENTS;
                case {'verbose', 'write_pdf', 'auto_contrast'}
                    validator = @islogical;
                case {'stipple_color','background_color'}
                    validator = @is_color_like;
                case 'white_cut'
                    validator = @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x <= 1;
                case 'supersample'
                    validator = @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1;
                case 'supersample_method'
                    validator = @(s) ischar(s) || isstring(s);
                case 'poisson_max_radius_factor'
                    validator = @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1;
                otherwise
                    error('Unknown parameter: %s', param_name);
            end
        end
        
        function is_valid = isValidImageFormat(filename)
            [~, ~, ext] = fileparts(filename);
            is_valid = any(strcmpi(ext, StipplerConfig.SUPPORTED_IMAGE_FORMATS));
        end
        
        function formats_str = getSupportedFormatsString()
            formats_str = strjoin(StipplerConfig.SUPPORTED_IMAGE_FORMATS, ', ');
        end
        
        function adaptive_eps = getAdaptiveEps(image_size)
            % Calculate adaptive epsilon based on image size
            % eps = BASE_EPS + (diagonal - 512√2) * EPS_SCALE_FACTOR / 100
            if nargin < 1 || isempty(image_size) || numel(image_size) ~= 2
                adaptive_eps = StipplerConfig.DEFAULT_EPS;
                return;
            end
            
            height = image_size(1);
            width = image_size(2);
            diagonal = sqrt(height^2 + width^2);
            reference_diagonal = sqrt(512^2 + 512^2);
            size_factor = (diagonal - reference_diagonal) / 100;
            adaptive_eps = StipplerConfig.BASE_EPS + size_factor * StipplerConfig.EPS_SCALE_FACTOR;
            adaptive_eps = max(StipplerConfig.MIN_EPS, min(adaptive_eps, StipplerConfig.MAX_EPS));
        end
        
        function enhanced_img = enhanceContrast(img, verbose)
            % Apply automatic contrast enhancement for low-contrast images
            if nargin < 2
                verbose = false;
            end
            
            img_std = std(img(:));
            
            if img_std < StipplerConfig.CONTRAST_THRESHOLD
                if verbose
                    fprintf('[INFO] Low contrast detected (std=%.3f), applying enhancement...\n', img_std);
                end
                
                img_min = min(img(:));
                img_max = max(img(:));
                img_range = img_max - img_min;
                
                if img_range > 1e-6
                    center = (img_min + img_max) / 2;
                    enhanced_img = center + (img - center) * StipplerConfig.CONTRAST_ENHANCEMENT_FACTOR;
                    enhanced_img = max(0, min(1, enhanced_img));
                    
                    if verbose
                        new_std = std(enhanced_img(:));
                        fprintf('[INFO] Contrast enhanced: std %.3f -> %.3f\n', img_std, new_std);
                    end
                else
                    enhanced_img = img;
                end
            else
                enhanced_img = img;
            end
        end
    end
end
