function [info, numberOfFrames,width,height] = getImageInfo(fname)

info = imfinfo(fname);
numberOfFrames = numel(info);
width = info.Width;
height = info.Height;

end
