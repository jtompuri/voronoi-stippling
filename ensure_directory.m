function ensure_directory(directory_path)
% ENSURE_DIRECTORY Create directory if it doesn't exist.
% directory_path: path to create; succeeds if already present.
% Throws an error on permission or other creation failures.

    if ~exist(directory_path, 'dir')
        try
            mkdir(directory_path);
        catch ME
            error('Failed to create directory "%s": %s', directory_path, ME.message);
        end
    end
end
