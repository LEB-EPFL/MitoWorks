function [FileName,PathName]=loadSegmentedImage()
    [FileName,PathName] = uigetfile('*.tif','Load the segmented image');
end
