classdef Stippler < handle
    % Core routines for weighted Voronoi stippling:
    % - Poisson-disk initialization (variable or fixed radius)
    % - Voronoi labeling via bwdist
    % - Weighted centroids via accumarray  
    % - Lloyd relaxation with early stop

    properties (Access = private)
        verbose
    end

    methods
        function obj = Stippler(verbose)
            if nargin < 1
                verbose = true;
            end
            obj.verbose = verbose;
            if obj.verbose
                fprintf('[INFO] MATLAB stippler initialized\n');
            end
        end

    function points = poisson_disk_sample(~, rho, n_points, min_dist_px, bias_density, maxRadiusFactor)
            % Poisson-disk sampling with variable or fixed radius
            if nargin < 5 || isempty(bias_density)
                bias_density = true;
            end
            if nargin < 6 || isempty(maxRadiusFactor)
                maxRadiusFactor = StipplerConfig.DEFAULT_POISSON_MAX_RADIUS_FACTOR;
            end
            [H,W] = size(rho);
            useVariable = (nargin < 4) || isempty(min_dist_px) || min_dist_px <= 0;
            mass = sum(rho(:));
            if useVariable
                r0 = max(1.0, sqrt(max(mass,1)/max(n_points,1)));
            else
                r0 = min_dist_px;
            end
            k = 30;
            cellSize = max(1, r0)/sqrt(2);
            gridW = ceil(W/cellSize);
            gridH = ceil(H/cellSize);
            grid = -ones(gridH, gridW, 'int32');
            pts = zeros(n_points*2, 2); nPts = 0;
            ptsR = zeros(n_points*2, 1); % local radius per point
            active = zeros(n_points*2, 1, 'int32'); nAct = 0;
            function [gi,gj] = gridIndex(x,y)
                gj = max(1, min(gridW, floor(x/cellSize)+1));
                gi = max(1, min(gridH, floor(y/cellSize)+1));
            end
            function rloc = localR(x,y)
                if useVariable
                    rr = rho(round(min(max(y,1),H)), round(min(max(x,1),W)));
                    rloc = r0 / sqrt(max(rr, 1e-6));
                    % Cap the local radius to avoid huge empty borders
                    rloc = min(rloc, r0 * maxRadiusFactor);
                else
                    rloc = r0;
                end
            end
            function v = reflectToRange(v, lo, hi)
                % Reflect a value into [lo, hi] (mirror at boundaries) to avoid
                % discarding candidates near borders. Handles multiple oversteps.
                if v >= lo && v <= hi
                    return;
                end
                period = 2*(hi - lo);
                if period <= 0
                    v = min(max(v, lo), hi);
                    return;
                end
                v = mod(v - lo, period);
                if v < 0, v = v + period; end
                if v > (hi - lo)
                    v = period - v;
                end
                v = v + lo;
            end
            function ok = farEnough(x,y,rc)
                [gi,gj] = gridIndex(x,y);
                % Expand neighbor range when rc is large
                range = 2 + ceil(rc/(r0+eps));
                for di=-range:range
                    ii = gi+di; if ii<1||ii>gridH, continue; end
                    for dj=-range:range
                        jj = gj+dj; if jj<1||jj>gridW, continue; end
                            idx = grid(ii,jj);
                            if idx>0
                                dx = pts(idx,1)-x; dy = pts(idx,2)-y;
                                % Required separation: max of local radii (stricter; prevents clumping at transitions)
                                req = max(rc, ptsR(idx));
                                if dx*dx+dy*dy < req*req, ok=false; return; end
                            end
                    end
                end
                ok = true;
            end
            function addPoint(x,y)
                nPts = nPts + 1;
                pts(nPts,:) = [x,y];
                ptsR(nPts) = localR(x,y);
                [gi,gj] = gridIndex(x,y);
                grid(gi,gj) = int32(nPts);
                nAct = nAct + 1; 
                active(nAct) = int32(nPts);
            end
            % Seed with a random valid point (avoid zero-density areas)
            max_tries = 10000; tries = 0;
            while tries < max_tries
                x0 = 1 + (W-1)*rand(); y0 = 1 + (H-1)*rand();
                rr0 = rho(round(y0), round(x0));
                % Never seed in strictly zero-density areas
                if rr0 <= 0
                    acc = false;
                elseif useVariable
                    % Optional biased seeding in variable mode when bias_density=true
                    acc = ~bias_density || (rand() < rr0);
                elseif ~useVariable && bias_density
                    % Fixed-radius: do not bias; uniform except zero-density
                    acc = true;
                else
                    acc = true;
                end
                if acc
                    addPoint(x0,y0); break;
                end
                tries = tries + 1;
            end
            % If no seed could be placed (e.g., zero-density image), return empty
            if nPts==0
                points = pts(1:0,:);
                return;
            end
            % Generate (with restarts)
            stuckRounds = 0;
            maxRelax = 12; % safety to avoid infinite relax loops
            while nPts<n_points
                if nAct==0
                    % restart from a new random seed
                    tries = 0; placed0 = false;
                    while tries < max_tries
                        x0 = 1 + (W-1)*rand(); y0 = 1 + (H-1)*rand();
                        rr0 = rho(round(y0), round(x0));
                        if rr0 <= 0
                            acc = false;
                        elseif useVariable
                            acc = true;
                        elseif ~useVariable && bias_density
                            acc = (rand() < rr0);
                        else
                            acc = true;
                        end
                        if acc
                            addPoint(x0,y0); placed0 = true; break;
                        end
                        tries = tries + 1;
                    end
                    if ~placed0
                        % Couldn't place a fresh seed: relax radius and retry
                        stuckRounds = stuckRounds + 1;
                        if stuckRounds <= maxRelax
                            % Shrink base radius to allow more placements
                            r0 = max(0.5, 0.95 * r0);
                            % Recompute per-point radii after r0 change
                            for qi=1:nPts
                                ptsR(qi) = localR(pts(qi,1), pts(qi,2));
                            end
                            % Rebuild acceleration grid at new cell size
                            cellSize = max(1, r0)/sqrt(2);
                            gridW = ceil(W/cellSize);
                            gridH = ceil(H/cellSize);
                            grid = -ones(gridH, gridW, 'int32');
                            for qi=1:nPts
                                [gi,gj] = gridIndex(pts(qi,1), pts(qi,2));
                                grid(gi,gj) = int32(qi);
                            end
                            % Explore with a bit less effort per active to avoid local traps
                            k = max(10, round(0.9 * k));
                            % Continue without breaking; try again to seed
                            continue;
                        else
                            break; % safety exit after too many relax attempts
                        end
                    end
                end
                % Random active with slight jitter to avoid repeatedly trying the same frontier
                ai = randi(nAct);
                pidx = active(ai);
                px = pts(pidx,1); py = pts(pidx,2);
                rA = ptsR(pidx);
                % For variable mode with clamp, take longer proposal steps based on unclamped radius
                if useVariable
                    rrA = rho(round(py), round(px));
                    rA_un = r0 / sqrt(max(rrA, 1e-6));
                    clampA = rA_un > (r0 * maxRadiusFactor);
                else
                    clampA = false; rA_un = rA;
                end
                placed = false;
                for t=1:k
                    % Propose in annulus; if clamped at active point, use unclamped step to explore further
                    if clampA
                        r = rA_un * (1+rand());
                    else
                        r = rA * (1+rand());
                    end
                    theta = 2*pi*rand();
                    x = px + r*cos(theta);
                    y = py + r*sin(theta);
                    % Reflect out-of-bounds candidates back into the domain to maintain
                    % consistent sampling efficiency near image borders.
                    x = reflectToRange(x, 1, W);
                    y = reflectToRange(y, 1, H);
                    if useVariable
                        rr = rho(round(y), round(x));
                        if rr <= 0, continue; end
                        % If destination area would be clamped (low density), add density gating when biasing
                        if bias_density
                            r_un = r0 / sqrt(max(rr,1e-6));
                            if r_un > (r0 * maxRadiusFactor) && rand() >= rr
                                continue;
                            end
                        end
                    else
                        % Fixed-radius: no density bias; only exclude zero-density
                        if rho(round(y), round(x)) <= 0, continue; end
                    end
                    rc = localR(x,y);
                    if farEnough(x,y,rc)
                        addPoint(x,y); placed = true; if nPts>=n_points, break; end
                    end
                end
                if ~placed
                    % remove from active by swapping with last
                    active(ai) = active(nAct); nAct = nAct - 1;
                end
            end
            points = pts(1:nPts,:);
        end

        function labels = compute_voronoi_labels(~, points, shape)
            % Voronoi labeling via bwdist
            height = shape(1);
            width = shape(2);

            if isempty(points)
                labels = zeros(height, width, 'uint32');
                return;
            end
            if ~exist('bwdist', 'file')
                error(['Image Processing Toolbox is required for Voronoi computation.\n' ...
                       'Please install the Image Processing Toolbox or use an alternative method.\n' ...
                       'Required function: bwdist']);
            end

            % Snap seeds to pixel centers and clamp
            xi = round(points(:, 1));
            yi = round(points(:, 2));
            xi = min(max(xi, 1), width);
            yi = min(max(yi, 1), height);

            % Seed label image: L(r,c) = seed index at that pixel, 0 elsewhere
            L = zeros(height, width, 'uint32');
            rc = sub2ind([height, width], yi, xi);

            % Handle duplicates by keeping first occurrence
            [rcu, ia] = unique(rc, 'stable');
            seed_ids = uint32((1:size(points, 1)).');
            L(rcu) = seed_ids(ia);

            % Nearest-seed labeling via distance transform
            BW = L > 0;
            [~, idx] = bwdist(BW, 'euclidean');
            labels = L(idx);
        end

        function new_points = compute_centroids(~, labels, rho, points)
            % Vectorized weighted centroids via accumarray
            [H, W] = size(labels);
            n_points = size(points, 1);

            lab = double(labels(:));
            w = double(rho(:));
            N = numel(lab);
            idx = (1:N)';
            % Correct mapping from linear index to (row, col) in MATLAB (column-major)
            rows = mod(idx - 1, H) + 1;           % y
            cols = floor((idx - 1) / H) + 1;      % x

            % Defensive clamp (should already be valid)
            cols = min(max(cols, 1), W);
            rows = min(max(rows, 1), H);

            sumw  = accumarray(lab, w,                 [n_points, 1], @sum, 0);
            sumwx = accumarray(lab, w .* double(cols), [n_points, 1], @sum, 0);
            sumwy = accumarray(lab, w .* double(rows), [n_points, 1], @sum, 0);

            new_points = points; % default to original where empty
            nz = sumw > 0;
            cx = zeros(n_points,1);
            cy = zeros(n_points,1);
            cx(nz) = sumwx(nz) ./ sumw(nz);
            cy(nz) = sumwy(nz) ./ sumw(nz);
            % Keep points at pixel-center integer coordinates
            new_points(nz, :) = [cx(nz), cy(nz)];
        end

    function [points, iters_done] = lloyd_relaxation(obj, points, rho, n_iter, verbose, eps)
            % Lloyd relaxation with early stopping based on mean displacement
            if nargin < 5
                verbose = obj.verbose;
            end
            if nargin < 6 || isempty(eps)
                eps = StipplerConfig.FALLBACK_EPS;
            end

            [height, width] = size(rho); %#ok<ASGLU>

            % Remove duplicate seeds before iterations
            [points, kept_idx] = obj.remove_duplicate_seeds(points, size(rho));
            if verbose && numel(kept_idx) < size(points,1)
                fprintf('[INFO] Removed %d duplicate seed(s). Now %d unique seeds.\n', ...
                    size(points,1) - numel(kept_idx), numel(kept_idx));
            end

            if verbose
                fprintf('[INFO] Lloyd relaxation with %d iterations (eps=%.3f) ...\n', n_iter, eps);
            end

            iters_done = 0;

            for iteration = 1:n_iter
                if verbose
                    fprintf('  → Iteration %d/%d\n', iteration, n_iter);
                end
                iter_start = tic;

                % Discrete Voronoi labeling and weighted centroid update
                labels = obj.compute_voronoi_labels(points, size(rho));
                new_points = obj.compute_centroids(labels, rho, points);

                % Calculate mean displacement for early stop
                disp_vec = hypot(new_points(:,1) - points(:,1), new_points(:,2) - points(:,2));
                mean_disp = mean(disp_vec);

                points = new_points;
                iters_done = iteration;

                if verbose
                    fprintf('    Completed in %.2f s (mean Δ = %.4f px)\n', toc(iter_start), mean_disp);
                end

                if mean_disp < eps
                    if verbose
                        fprintf('    Early stop: mean displacement < eps (%.4f < %.4f)\n', mean_disp, eps);
                    end
                    break;
                end
            end

            if verbose
                fprintf('[INFO] Lloyd relaxation complete after %d iterations\n', iters_done);
            end
        end

        function [points_out, kept_idx] = remove_duplicate_seeds(~, points, shape)
            % Remove duplicate seeds that snap to the same pixel after rounding.
            % Keeps first occurrence (stable), returns filtered points and indices kept.
            height = shape(1);
            width  = shape(2);
            xi = round(points(:, 1));
            yi = round(points(:, 2));
            xi = min(max(xi, 1), width);
            yi = min(max(yi, 1), height);
            rc = sub2ind([height, width], yi, xi);
            [~, ia] = unique(rc, 'stable');
            points_out = points(ia, :);
            kept_idx = ia;
        end

    end
end
