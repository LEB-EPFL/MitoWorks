function [PullingFrames1, PullingFrames2]=splitMito(BreakingFrame,BreakFrame,EndFrame,SegmentedMitoIDImageDM1,SegmentedMitoIDImageDM2);
% function divides the segmented image in the frame before fission and
% allows the tracking of both daughter mitochondria after fission
% 
%   INPUTS:
%   BreakingFrame = last frame before fission
%   EndFrame = last frame
%   SegmentedMitoIDImageDM1 and DM2 = segmented images of daughter
%   mitochondria after fission
% 
%   OUTPUTS:
%   PullingFrames1 and 2 = segmented images including the last frame before
%   fission where mitochondria are "artificially" broken

    % segment and label image
    ccp = bwconncomp(BreakingFrame, 4);
    
    % isolate daughter mito
    grain1 = false(size(BreakingFrame));
    grain1(ccp.PixelIdxList{1}) = true;
    grain2 = false(size(BreakingFrame));
    grain2(ccp.PixelIdxList{2}) = true;
    
    % assign to right frames
    PullingFrames1{BreakFrame}=grain1;
    PullingFrames2{BreakFrame}=grain2;

    for f=(BreakFrame+1):EndFrame
        PullingFrames1{f}=SegmentedMitoIDImageDM1{f};
        PullingFrames2{f}=SegmentedMitoIDImageDM2{f};
    end
    
end
