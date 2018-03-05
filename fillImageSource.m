function [pixel_size, limit_speed, filter_area] = fillImageSource(image_source)

% pixel_size = size of pixel [nm]
% limit_distance = filter maximum distance between frames for same mito [nm]
% filter_area = minimum area of mito [pixels], usually around 50000nm^2 

if strcmp(image_source,'SIM')==1
    pixel_size=30; limit_speed=200; filter_area=50000/900;
elseif strcmp(image_source,'ZEISS')==1
    pixel_size=100; limit_speed=200; filter_area=5;
end

end

