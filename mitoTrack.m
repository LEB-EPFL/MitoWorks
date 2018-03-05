function [Speeds, Tracks, SegmentedMitoIDImage] = mitoTrack(Fname,Timestep,IndividualMito,IndividualMitoID,StartFrame,EndFrame,PixelSize)
% Tracks centres of mass of mitochondria to output their speeds [nm],
% tracks [in pixels] and a segmented image of a single mitochondrion if
% IndiviualMito = true
%
% Inputs:
% Fname = Path to segmented tif image
% Timestep = time between frames [s]
% IndividualMito = true if you would like to just follow 1 mito
% IndividualMitoID = ID of track which to follow
% StartFrame = number of first frame to analyse
% EndFrame = number of last frame to analyse
% PixelSize = pixel size [nm]
% 


% [PathStr,Name,Ext] = fileparts(Fname);


ImageSource='SIM';

LimitSpeed= (4000/PixelSize)/Timestep;
LimitDistance=LimitSpeed*Timestep;

    h = waitbar(0,'Running mitoTrack.m');

nTracks=0;
ImageBW=cell(1,(EndFrame));
CC=cell(1,(EndFrame));
LabeledImageBW=cell(1,(EndFrame));
MitoMeasurements=cell(1,(EndFrame));
NumberOfMito=zeros(1,(EndFrame));
Centroids=cell(1,(EndFrame));
IDX=cell(1,(EndFrame));
DistancesMito=cell(1,(EndFrame));
nonAssignedIdx=cell(1,(EndFrame));
% Speeds
for f=StartFrame:EndFrame
    waitbar((f-StartFrame)/(EndFrame-StartFrame));
    ImageBW{f} = convertBinaryImage(Fname,f);
    [CC{f}, LabeledImageBW{f}] = labelBWImage(ImageBW{f},4);
    MitoMeasurements{f} = regionprops(LabeledImageBW{f}, ImageBW{f}, 'all');
    NumberOfMito(f)=size(MitoMeasurements{f},1);
    
    % isolate centroids`
    for i=1:NumberOfMito(f)
        Centroids{f}(i,:)=MitoMeasurements{f}(i).Centroid;
        
        % initialize tracks
        if f==StartFrame
            nTracks=nTracks+1;
            MitoMeasurements{f}(i).Track=nTracks;
            MitoMeasurements{f}(i).StartFrame=StartFrame;
            Tracks{nTracks}(StartFrame,:)=MitoMeasurements{f}(i).Centroid;
            StartFrames{nTracks}=MitoMeasurements{f}(i).StartFrame;
        end
    end
    
    if f>StartFrame
        CentroidsBefore=Centroids{f-1};
        CentroidsNow=Centroids{f};       
        
        % find mito in next frame that plot to ones in previous frame
        
        % IDX = knnsearch(X,Y) finds the nearest neighbor in X for each 
        % point in Y. IDX is a column vector with my rows. Each row in IDX 
        % contains the index of nearest neighbor in X for the corresponding 
        % row in Y.
        [IDX{f},DistancesMito{f}]=knnsearch(CentroidsNow,CentroidsBefore);
        
        idx=IDX{f};
        solvedIDX=zeros(1,max(idx));

        
        for i=1:numel(IDX{f})
        
            if((numel(idx(idx==idx(i)))>1) && solvedIDX(idx(i))==0)
                % find shortest distance for multiple matches
                distances=DistancesMito{f};
                distances=distances(find(idx==idx(i)));
                minD=min(distances);
                distances=DistancesMito{f};
                minidx=find(idx==idx(i) & distances==minD);
%                 minidx=find(DistancesMito{f}(find(idx==idx(i)))==min(distances));
                if (minD<LimitDistance)
                    otheridx=setdiff(find(idx==idx(i)),minidx);
                    MitoMeasurements{f}(idx(minidx)).Track=MitoMeasurements{f-1}(minidx).Track;
                    MitoMeasurements{f}(idx(minidx)).Speed=sqrt((MitoMeasurements{f}(idx(minidx)).Centroid(1)- MitoMeasurements{f-1}(minidx).Centroid(1))^2 + (MitoMeasurements{f}(idx(minidx)).Centroid(2) - MitoMeasurements{f-1}(minidx).Centroid(2))^2)*PixelSize/Timestep; 
                    MitoMeasurements{f}(idx(minidx)).StartFrame=MitoMeasurements{f-1}(minidx).StartFrame;

                    Speeds(f,MitoMeasurements{f}(idx(minidx)).Track)=MitoMeasurements{f}(idx(minidx)).Speed;
                    StartFrames{MitoMeasurements{f}(idx(minidx)).Track}=MitoMeasurements{f}(idx(minidx)).StartFrame;
                    Tracks{MitoMeasurements{f}(idx(minidx)).Track}(f,:)=MitoMeasurements{f}(idx(minidx)).Centroid;
%                     
                else
                    otheridx=find(idx==idx(i));
                end
                    
                solvedIDX(idx(i))=1;
%                 for j=1:numel(otheridx)
%                     nTracks=nTracks+1;
%                     MitoMeasurements{f}(idx(otheridx(j))).Track=nTracks;
%                     MitoMeasurements{f}(idx(otheridx(j))).StartFrame=f;
%                     StartFrames{MitoMeasurements{f}(idx(otheridx(j))).Track}=f;
%                     Tracks{MitoMeasurements{f}(idx(otheridx(j))).Track}(f,:)=MitoMeasurements{f}(idx(otheridx(j))).Centroid;
%                 end
%                 
            elseif(solvedIDX(idx(i))==0)
                % assign matched mito to track
                if DistancesMito{f}(i)<LimitDistance
                    MitoMeasurements{f}(idx(i)).Track=MitoMeasurements{f-1}(i).Track;
                    MitoMeasurements{f}(idx(i)).Speed=sqrt((MitoMeasurements{f}(idx(i)).Centroid(1)- MitoMeasurements{f-1}(i).Centroid(1))^2 + (MitoMeasurements{f}(idx(i)).Centroid(2) - MitoMeasurements{f-1}(i).Centroid(2))^2)*PixelSize/Timestep; 
                    Speeds(f,MitoMeasurements{f-1}(i).Track)=MitoMeasurements{f}(idx(i)).Speed;
                    MitoMeasurements{f}(idx(i)).StartFrame=MitoMeasurements{f-1}(i).StartFrame;
                    Tracks{MitoMeasurements{f}(idx(i)).Track}(f,:)=MitoMeasurements{f}(idx(i)).Centroid;
                    StartFrames{MitoMeasurements{f}(idx(i)).Track}=MitoMeasurements{f}(idx(i)).StartFrame;
                else
                    nTracks=nTracks+1;
                    MitoMeasurements{f}(idx(i)).Track=nTracks;
                    MitoMeasurements{f}(idx(i)).StartFrame=f;
                    Tracks{MitoMeasurements{f}(idx(i)).Track}(f,:)=MitoMeasurements{f}(idx(i)).Centroid;
                    StartFrames{MitoMeasurements{f}(idx(i)).Track}=f;
                end
            end
            
        end
        
        % assign unassigned mito to new tracks
        nonAssignedIdx{f}=setdiff([1:NumberOfMito(f)],IDX{f});
        
        for i=1:numel(nonAssignedIdx{f})
            nTracks=nTracks+1;
            MitoMeasurements{f}(nonAssignedIdx{f}(i)).Track=nTracks;
            MitoMeasurements{f}(nonAssignedIdx{f}(i)).StartFrame=f;
            Tracks{MitoMeasurements{f}(nonAssignedIdx{f}(i)).Track}(f,:)=MitoMeasurements{f}(nonAssignedIdx{f}(i)).Centroid;
            StartFrames{MitoMeasurements{f}(nonAssignedIdx{f}(i)).Track}=MitoMeasurements{f}(nonAssignedIdx{f}(i)).StartFrame;
        end
    end
    
    idxIM(f)=0;
    for i=1:NumberOfMito(f)
        if MitoMeasurements{f}(i).Track==IndividualMitoID & IndividualMito==true
            idxIM(f)=i;
        end
    end
    
end
close(h);

EndFrames=cell(1,length(nTracks));
for j=1:nTracks
    EndFrames{j}=size(Tracks{j},1);
end

clear Centroids;

SegmentedMitoIDImage=cell(1,(EndFrame));
for f=StartFrame:EndFrame
    FrameTracks=0;
    for i=1:NumberOfMito(f)
        for j=1:nTracks
            if StartFrames{j}<=f && EndFrames{j}>=f
                FrameTracks=FrameTracks+1;
                Centroids{f}(FrameTracks,1:2)=Tracks{j}(f,:);
                Centroids{f}(FrameTracks,3)=j;
            end
        end
    end
    if IndividualMito==true
        SpecificLabel=false(size(ImageBW{f}));
        if idxIM(f) >0
            SpecificLabel(CC{f}.PixelIdxList{idxIM(f)})=true;
        end
%         figure
        imshow(SpecificLabel); 
        SegmentedMitoIDImage{f}=SpecificLabel;
        hold on
            if StartFrames{IndividualMitoID}<=f && EndFrames{IndividualMitoID}>=f
                line(Tracks{IndividualMitoID}(StartFrames{IndividualMitoID}:EndFrames{IndividualMitoID},1),Tracks{IndividualMitoID}(StartFrames{IndividualMitoID}:EndFrames{IndividualMitoID},2),'LineWidth', 2);
                plot(Tracks{IndividualMitoID}(f,1),Tracks{IndividualMitoID}(f,2),'ro');
            end
        hold off
        title(sprintf('Frame %d',f));drawnow;
%         averageSpeed=mean(nonzeros(Speeds(:,IndividualMitoID)));
%         stdSpeed=std(nonzeros(Speeds(:,IndividualMitoID)));

    else
%         figure
        labeled=labelmatrix(CC{f});
        RBG_label=label2rgb(labeled,'hsv','k','shuffle');
        hold on

        RGB=insertText(RBG_label,Centroids{f}(:,1:2),Centroids{f}(:,3),'AnchorPoint','Center','FontSize',20,'BoxOpacity',0.8);
        imshow(RGB);
    	SegmentedMitoIDImage{f}=RGB;

        hold off
        for j=1:nTracks
            if StartFrames{j}<=f && EndFrames{j}>=f
                hold on
                line(Tracks{j}(StartFrames{j}:EndFrames{j},1),Tracks{j}(StartFrames{j}:EndFrames{j},2),'LineWidth', 2);
                plot(Tracks{j}(f,1),Tracks{j}(f,2),'ro');
                hold off
            end
        end
        title(sprintf('Frame %d',f));drawnow;
        averageSpeed=mean(nonzeros(Speeds));
        stdSpeed=std(nonzeros(Speeds));
    end
end



% figure
if IndividualMito==true
    plot(find(Speeds(:,IndividualMitoID)),nonzeros(Speeds(:,IndividualMitoID)));
    xlabel('Frame number');
    ylabel('Speed [nm/s]');
else

    for t=1:size(Speeds,2)
        hold on
        plot(find(Speeds(:,t)),nonzeros(Speeds(:,t)));
        hold off
        
    end
    xlabel('Frame number');
    ylabel('Speed [nm/s]');
end
end
