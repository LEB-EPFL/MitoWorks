function ComposedSegImage = combineParentDaughterSegImage(Orig,DM1,DM2,StartFrame,EndFrame,BreakFrame);

    for f=StartFrame:BreakFrame-1
        ComposedSegImage{f}=Orig{f};
    end

    for f=BreakFrame:EndFrame
        ComposedSegImage{f}=DM1{f}+DM2{f};
    end
end
