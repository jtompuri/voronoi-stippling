function tf = is_color_like(c)
% IS_COLOR_LIKE True if input is a MATLAB color spec (name or [r g b])
% Accepts char/string names or numeric 1x3 RGB in [0,1]

    if ischar(c) || isstring(c)
        tf = true; % defer name validation to graphics layer
        return;
    end
    tf = isnumeric(c) && numel(c) == 3 && all(isfinite(c(:))) && all(c(:) >= 0) && all(c(:) <= 1);
end
