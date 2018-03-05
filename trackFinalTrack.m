function [SelectedTracksSpeed]=trackFinalTrack(DrpTrackCons,TimeScale)

if size(DrpTrackCons,1)>1
    for f=2:size(DrpTrackCons,1)
        SelectedTracksSpeed(1,DrpTrackCons(f,4))=sqrt((DrpTrackCons(f,2)-DrpTrackCons(f-1,2))^2+(DrpTrackCons(f,3)-DrpTrackCons(f-1,3))^2)/TimeScale;
    end
else
    SelectedTracksSpeed(1,DrpTrackCons(1,4))=0;
end

