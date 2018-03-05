function [Drp1Tracks, RadiusTrack ] = FilterTracks(Drp1TrackFile, StartFrame, EndFrame)
Tracks=importdata(Drp1TrackFile);
TracksImported=Tracks.data;
ColNames=Tracks.textdata;

[~, Col_TrackID]=find(strcmp('TRACK_ID', ColNames));
[~, Col_X]=find(strcmp('POSITION_X', ColNames));
[~, Col_Y]=find(strcmp('POSITION_Y', ColNames));
% [row_Time Col_Time]=find(strcmp('POSITION_T', ColNames));
[~, Col_Frame]=find(strcmp('FRAME', ColNames));
[~, Col_R]=find(strcmp('RADIUS', ColNames));


%text and data are shifted by 2 columns
Drp1Tr(:,1)=TracksImported(:,(Col_TrackID-2))+1;
Drp1Tr(:,2)=TracksImported(:,(Col_X-2) );
Drp1Tr(:,3)=TracksImported(:,(Col_Y-2) );
Drp1Tr(:,4)=TracksImported(:,(Col_Frame-2))+1;
Drp1Tr(:,5)=TracksImported(:,(Col_R-2));

RadiusTrack=Drp1Tr(1,5);


%separate each track ID
count=1;
Drp1Tracks_all=cell(max(Drp1Tr(:,1)),max(Drp1Tr(:,1)));
for i=1:max(Drp1Tr(:,1));
target=find(Drp1Tr(:,1)==i);
Drp1Tracks_all{count,1}(:,1)=Drp1Tr(target,1);%track ID
Drp1Tracks_all{count,1}(:,2)=Drp1Tr(target,2);%x position
Drp1Tracks_all{count,1}(:,3)=Drp1Tr(target,3);%y position
Drp1Tracks_all{count,1}(:,4)=Drp1Tr(target, 4);%frame
count=count+1;
end

%remove empty cell contents
Drp1Tracks_all=Drp1Tracks_all(~cellfun('isempty', Drp1Tracks_all)); 

%find tracks within specified frame range
ValidTrackIDrow=cell(1,length(Drp1Tracks_all));
ValidTrackIDcol=cell(1,length(Drp1Tracks_all));
Drp1Tracks=cell(length(Drp1Tracks_all),1);
for t=1:length(Drp1Tracks_all)
    
   [ValidTrackIDrow{t}, ValidTrackIDcol{t}]=find(Drp1Tracks_all{t,1}(:,4)>=StartFrame & Drp1Tracks_all{t,1}(:,4)<=EndFrame);
   Drp1Tracks{t,1}=Drp1Tracks_all{t,1}( ValidTrackIDrow{1,t}, :);
end

Drp1Tracks=Drp1Tracks(~cellfun('isempty', Drp1Tracks)); 

end

