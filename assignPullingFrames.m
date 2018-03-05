function PullingFrames = assignPullingFrames(BreakingFrame,SegmentedMitoIDImageDM1,SegmentedMitoIDImageDM2,BreakFrame,EndFrame);
PullingFrames=cell(1,EndFrame);
    PullingFrames{BreakFrame}=BreakingFrame;
    
    for f=(BreakFrame+1):EndFrame
        PullingFrames{f}=max(SegmentedMitoIDImageDM1{f},SegmentedMitoIDImageDM2{f});
    end

end
