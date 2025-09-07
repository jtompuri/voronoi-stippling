% filepath: save_as_tsp.m
function save_as_tsp(points, filename, name)
% SAVE_AS_TSP Save stipples as a TSPLIB TSP file.
%
% Parameters:
%   points    Nx2 [x,y] coordinates in image pixels
%   filename  output path (.tsp)
%   name      TSPLIB NAME field

    % Input validation
    if ~isnumeric(points) || size(points, 2) ~= 2
        error('points must be an Nx2 numeric array');
    end
    if ~ischar(filename) && ~isstring(filename)
        error('filename must be a string or character array');
    end
    if ~ischar(name) && ~isstring(name)
        error('name must be a string or character array');
    end
    if size(points, 1) == 0
        warning('No points to save in TSP file');
    end

    % Ensure output directory exists
    [output_dir, ~, ~] = fileparts(filename);
    ensure_directory(output_dir);
    
    fprintf('[INFO] Saving stipples to TSP: %s\n', filename);
    
    % Open file for writing
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot write TSP file: %s', filename);
    end
    
    try
        % Write TSP header
        fprintf(fid, 'NAME: %s\n', name);
        fprintf(fid, 'TYPE: TSP\n');
        fprintf(fid, 'COMMENT: Stipple points for Traveling Salesman Problem\n');
        fprintf(fid, 'DIMENSION: %d\n', size(points, 1));
        fprintf(fid, 'EDGE_WEIGHT_TYPE: EUC_2D\n');
        fprintf(fid, 'NODE_COORD_SECTION\n');
        
        % Write coordinates (1-indexed nodes)
        for i = 1:size(points, 1)
            fprintf(fid, '%d %.6f %.6f\n', i, points(i, 1), points(i, 2));
        end
        
        fprintf(fid, 'EOF\n');
        fclose(fid);
        
    catch ME
        fclose(fid);
        error('Error writing TSP file: %s', ME.message);
    end
end
