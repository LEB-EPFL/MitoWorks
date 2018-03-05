function [SelectedTracksSpeed] = trackSelectedTracks(SelectedTracks,Timescale)

for t=1:size(SelectedTracks,1)
    if size(SelectedTracks{t},1)>1
        for f=2:size(SelectedTracks{t},1)
            SelectedTracksSpeed(t,SelectedTracks{t}(f,4))=sqrt((SelectedTracks{t}(f,2)-SelectedTracks{t}(f-1,2))^2+(SelectedTracks{t}(f,3)-SelectedTracks{t}(f-1,3))^2)/Timescale;
        end
    else
        SelectedTracksSpeed(t,SelectedTracks{t}(1,4))=0;
    end
end