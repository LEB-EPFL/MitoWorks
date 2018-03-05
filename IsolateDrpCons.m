function [ DrpTrackCons] = IsolateDrpCons( SelectedTracks,  Drp1EndFrame, IDmin, subNormals, Pixel )

%output is the Drp1 track at the constriction diameter, which leads to
%fission

   %find location of min diameter on image itself
   IDminDrp1=IDmin(Drp1EndFrame);
   X_minD=mean(subNormals{1, Drp1EndFrame}(IDminDrp1,1:2:3));
   Y_minD=mean(subNormals{1, Drp1EndFrame}(IDminDrp1,2:2:4));
   MinPoint=[X_minD, Y_minD];
   
   %find candidates for track of interest, i.e those tracks with a
   %localization at the final frame
   SpotsOfInterest=cell(size(SelectedTracks, 1),1);
   for i =1:size(SelectedTracks, 1)
       for j=1:size(SelectedTracks{i,1}, 1)
       if SelectedTracks{i, 1}(j,4)==Drp1EndFrame;
       SpotsOfInterest{i,1}(j, :)=SelectedTracks{i,1}(j,:);  
       end
       end
   end
   SpotsOfInterest=cell2mat(SpotsOfInterest);
      SpotsOfInterest(~any(SpotsOfInterest,2), :)=[]; %delete rows with only zeros
   
   %find the spot of interest closest to the min diameter
idx=knnsearch(SpotsOfInterest(:,2:1:3)./Pixel, MinPoint)

TrackIDOfInterest=SpotsOfInterest(idx, 1); %ID of the track at the constriction site

%insert entire tracks info into a cell
DrpTrackConsID=cell(size(SelectedTracks, 1),1);
for m =1:size(SelectedTracks, 1)
       
           if SelectedTracks{m, 1}(1,1)==TrackIDOfInterest;
               DrpTrackConsID{m,1}=SelectedTracks{m,1};
           end

end
 DrpTrackCons=cell2mat(DrpTrackConsID);
