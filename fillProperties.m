function [visualisation, limit_shape_change, with_area] = fillProperties

% visualisation = CM/BB visualisation during processing (false faster)
% limit_shape_change = filter maximum possible shape change [percentage]
% with_area = consider area for mito tracking

visualisation = true;

limit_shape_change = 70;

with_area = false;

end