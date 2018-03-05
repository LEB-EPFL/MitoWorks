function [frameStamp_multiTracks] = frameStamper(SelectedTracks, FrameBeforeFission)


    for i=1:size(SelectedTracks,1)
        frames=zeros(size(SelectedTracks{i,1}(:,4),1),1);
        for j=1: size(SelectedTracks{i,1}(:,4),1)
            if i==1
              frames(j,1)=SelectedTracks{i,1}(j,4);   
            end
            if i>1
%     frames(j+size(SelectedTracks{i-1,1}(:,4),1),1)=SelectedTracks{i,1}(j,4);
frames(j+size(frames,1))=SelectedTracks{i,1}(j,4);
            end
        end
    end
    
 
    OrderFrames=sort(frames,'ascend');
    frameStamp_multiTracks=unique(OrderFrames);
    
    for p=1:size(frameStamp_multiTracks,1)
if    frameStamp_multiTracks(p)> FrameBeforeFission
    frameStamp_multiTracks(p)=NaN;
    
end
if    frameStamp_multiTracks(p)==0
    frameStamp_multiTracks(p)=NaN;
end

    end
     frameStamp_multiTracks= frameStamp_multiTracks(~isnan(frameStamp_multiTracks));
end

