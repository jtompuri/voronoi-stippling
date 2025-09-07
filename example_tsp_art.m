% example_tsp_art.m - TSP tour visualization driver

% Configuration
exact_basename = 'example-1024px_10000';
base_prefix    = '';

% Curve options
use_spline     = false;
spline_method  = 'pchip';
spline_samples = [];

cfg = StipplerConfig.getDefault();
cfg.stipple_color = 'black';
cfg.background_color = 'white';

tsp_dir  = fullfile('stipplings','tsp');
tour_dir = fullfile('stipplings','tours');
if ~exist(tour_dir,'dir')
    tour_dir = tsp_dir;
end

pick = '';
function tf = has_tour(basename, tour_dir)
    tf = exist(fullfile(tour_dir, [basename '.tour']), 'file') ~= 0 || ...
         exist(fullfile(tour_dir, [basename '.lk']),   'file') ~= 0 || ...
         exist(fullfile(tour_dir, [basename '.sol']),  'file') ~= 0;
end

if ~isempty(exact_basename)
    tsp_exact = fullfile(tsp_dir, [exact_basename '.tsp']);
    if exist(tsp_exact,'file')
        pick = exact_basename;
    else
        warning('Exact basename not found: %s', tsp_exact);
    end
end

if isempty(pick)
    listing = dir(fullfile(tsp_dir, sprintf('%s_*.tsp', base_prefix)));
    if isempty(listing)
        error('No TSP files found for prefix: %s in %s', base_prefix, tsp_dir);
    end
    candidates = struct('basename',{},'num',{},'datenum',{},'hasTour',{});
    for k = 1:numel(listing)
        name = listing(k).name;
        [~, bn, ~] = fileparts(name);
        tok = regexp(bn, sprintf('^%s_(\\d+)$', regexptranslate('escape', base_prefix)), 'tokens', 'once');
        if isempty(tok), continue; end
        num = str2double(tok{1});
        candidates(end+1) = struct( ...
            'basename', bn, ...
            'num', num, ...
            'datenum', listing(k).datenum, ...
            'hasTour', has_tour(bn, tour_dir)); %#ok<SAGROW>
    end
    if isempty(candidates)
        error('No candidates with numeric suffix found for prefix: %s', base_prefix);
    end
    withTour = candidates([candidates.hasTour]);
    if ~isempty(withTour)
        [~, idx] = max([withTour.num]);
        pick = withTour(idx).basename;
    else
        [~, idx] = max([candidates.num]);
        pick = candidates(idx).basename;
    end
end

fprintf('[INFO] Visualizing basename: %s\n', pick);

% Validate and find matching TSP/tour pair
[tsp_ok, detail] = validate_pair(pick, tsp_dir, tour_dir);
if ~tsp_ok
    warning('[WARN] %s', detail);
    alt = find_matching_pair(base_prefix, tsp_dir, tour_dir);
    if ~isempty(alt)
        fprintf('[INFO] Switching to candidate with matching counts: %s\n', alt);
        pick = alt;
    else
        error('No matching TSP/tour pair found for prefix %s', base_prefix);
    end
end

sel_tsp  = fullfile(tsp_dir,  [pick '.tsp']);
sel_tour = first_existing_tour(pick, tour_dir);
if isempty(sel_tour)
    error('No tour file found for %s in %s', pick, tour_dir);
end
fprintf('[INFO] Using TSP:  %s\n', sel_tsp);
fprintf('[INFO] Using TOUR: %s\n', sel_tour);

visualize_tsp_tour(pick, ...
    'tsp_path', sel_tsp, ...
    'tour_path', sel_tour, ...
    'line_width', 2.0, ...
    'line_color', cfg.stipple_color, ...
    'background_color', cfg.background_color, ...
    'join_style', 'miter', ...
    'use_spline', use_spline, ...
    'spline_method', spline_method, ...
    'spline_samples', spline_samples);

% =========================
% Helper functions
% =========================
function [ok, detail] = validate_pair(basename, tsp_dir, tour_dir)
    tsp_path  = fullfile(tsp_dir,  [basename '.tsp']);
    tour_path = first_existing_tour(basename, tour_dir);
    if ~exist(tsp_path,'file')
        ok = false; detail = sprintf('Missing TSP: %s', tsp_path); return; end
    if isempty(tour_path) || ~exist(tour_path,'file')
        ok = false; detail = sprintf('Missing tour for basename: %s', basename); return; end
    ncoords = read_tsp_coord_count(tsp_path);
    [hdr_count, max_idx] = read_tour_header_and_max(tour_path);
    if isfinite(hdr_count) && hdr_count > 0 && hdr_count ~= ncoords
        ok = false; detail = sprintf('Header count %d \x2260 TSP coords %d', hdr_count, ncoords); return; end
    if max_idx + 1 > ncoords
        ok = false; detail = sprintf('Tour max index %d exceeds TSP coords %d', max_idx, ncoords); return; end
    ok = true; detail = 'OK';
end

function best = find_matching_pair(prefix, tsp_dir, tour_dir)
    best = '';
    listing = dir(fullfile(tsp_dir, sprintf('%s_*.tsp', prefix)));
    if isempty(listing), return; end
    candidates = struct('bn',{},'ncoords',{},'hdr',{},'maxidx',{});
    for k = 1:numel(listing)
        [~, bn, ~] = fileparts(listing(k).name);
        tsp_path  = fullfile(tsp_dir, [bn '.tsp']);
        tour_path = first_existing_tour(bn, tour_dir);
        if ~exist(tsp_path,'file') || isempty(tour_path) || ~exist(tour_path,'file')
            continue;
        end
        ncoords = read_tsp_coord_count(tsp_path);
        [hdr_count, max_idx] = read_tour_header_and_max(tour_path);
        candidates(end+1) = struct('bn',bn,'ncoords',ncoords,'hdr',hdr_count,'maxidx',max_idx); %#ok<SAGROW>
    end
    if isempty(candidates), return; end
    % Match rule priority: (1) hdr == ncoords; (2) maxidx+1 == ncoords; (3) largest ncoords
    exact = find([candidates.hdr] == [candidates.ncoords], 1, 'first');
    if ~isempty(exact)
        best = candidates(exact).bn; return; end
    idx2 = find(([candidates.maxidx] + 1) == [candidates.ncoords], 1, 'first');
    if ~isempty(idx2)
        best = candidates(idx2).bn; return; end
    % Fallback to the largest coordinate count (likely the newest filtered run)
    [~, idx3] = max([candidates.ncoords]);
    best = candidates(idx3).bn;
end

function tp = first_existing_tour(basename, tour_dir)
    exts = {'.tour','.lk','.sol'};
    tp = '';
    for i = 1:numel(exts)
        cand = fullfile(tour_dir, [basename exts{i}]);
        if exist(cand,'file'), tp = cand; return; end
    end
end

function n = read_tsp_coord_count(tsp_path)
    n = 0;
    fid = fopen(tsp_path, 'r'); if fid == -1, return; end
    c = onCleanup(@() fclose(fid)); 
    while ~feof(fid)
        line = fgetl(fid); if ~ischar(line), break; end
        if contains(lower(line), 'node_coord_section')
            % Count the lines until EOF or another section marker
            while ~feof(fid)
                pos = ftell(fid);
                l2 = fgetl(fid);
                if ~ischar(l2) || isempty(l2) || contains(lower(l2),'eof') || contains(lower(l2),'display_data_section')
                    break;
                end
                vals = sscanf(l2, '%d %f %f');
                if numel(vals) == 3
                    n = n + 1;
                else
                    fseek(fid, pos, 'bof'); break;
                end
            end
            break;
        end
    end
end

function [hdr, mx] = read_tour_header_and_max(tour_path)
    hdr = NaN; mx = -Inf;
    fid = fopen(tour_path, 'r'); if fid == -1, return; end
    c = onCleanup(@() fclose(fid)); 
    line = fgetl(fid);
    if ischar(line)
        hc = str2double(strtrim(line));
        if isfinite(hc) && hc > 0
            hdr = hc;
        end
    end
    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        nums = sscanf(strtrim(line), '%d %d');
        if numel(nums) >= 2
            mx = max(mx, max(nums));
        end
    end
    if ~isfinite(mx), mx = -1; end
end
