function [ SelectedTracks] = SelectTracks(MitoImage,DrpImage, StartFrame, EndFrame, Drp1Tracks, Pixel )

figure
imshow(imfuse(imread(MitoImage, StartFrame),(imread(DrpImage, StartFrame))), [])
hold on
%plot an asteriks at the start of the track
c=linspace(1,10*length(Drp1Tracks),length(Drp1Tracks));
   area=30;
   h=waitbar(10, 'Please wait. Plotting tracks');
for i=1:length(Drp1Tracks)
hold on
title('All tracks plotted on start frame. Draw a polygon around tracks (solid dots) of interest. ')
colormap(jet)
    scatter(Drp1Tracks{i, 1}(1,2)./Pixel,Drp1Tracks{i, 1}(1,3)./Pixel,area,c(:,i),'filled') 
    for j=1:length(Drp1Tracks{i,1})
        plot(Drp1Tracks{i, 1}(:,2)./Pixel,Drp1Tracks{i, 1}(:,3)./Pixel, 'c-')
    end
    waitbar(i/length(Drp1Tracks));
   
end
close(h);
%user draws polygon around tracks
    hPoly=impoly;
   PosPoly=hPoly.getPosition();
   close
   
%find start of tracks inside user-defined polygon 
inPoly=zeros(length(Drp1Tracks),1);
onPoly=zeros(length(Drp1Tracks),1);
SelectTracks=cell(length(Drp1Tracks),1);
SelectedTracks=cell(length(Drp1Tracks),1);
for l=1:length(Drp1Tracks)
   [inPoly(l,1), onPoly(l,1)]=inpolygon(Drp1Tracks{l,1}(1,2)./Pixel, Drp1Tracks{l,1}(1,3)./Pixel, PosPoly(:,1), PosPoly(:,2));
%    [TOIr{l,1}, TOIc{l,1}]=find(inPoly(l,1)==1|onPoly(l,1)==1);
    TOI=inPoly+onPoly;
    TOI=(TOI>0);
   SelectTracks{l,1}=TOI(l,1)*Drp1Tracks{l,1};
   
   if TOI(l,1)>0
   SelectedTracks{l,1}=SelectTracks{l,1};
   end
   
end
 %remove empty cells
   SelectedTracks=SelectedTracks(~cellfun('isempty', SelectedTracks)); 

    %plot the track (dot localization) for each frame with label ID
    c2=linspace(1,10*length(SelectedTracks),length(SelectedTracks));
for f=StartFrame:EndFrame
       figure;
       imshow(imfuse(imread(MitoImage, f),(imread(DrpImage, f))), []);
        str1 = sprintf('frame %d ', f);
    title(str1)
       hold on
       for n=1:length(SelectedTracks)
    
           for m=1:size(SelectedTracks{n,1}, 1)
           if SelectedTracks{n,1}(m,4)==f
        scatter(SelectedTracks{n, 1}(m, 2)./Pixel,SelectedTracks{n, 1}(m, 3)./Pixel, area,c2(n), 'filled');
        str2=sprintf('ID %d', SelectedTracks{n, 1}(m, 1));
        text(SelectedTracks{n, 1}(m, 2)./Pixel+3,SelectedTracks{n, 1}(m, 3)./Pixel+3, str2);
           end
           end
       end
end

