function [TracksinFrames I] = measureDrpSignal_multiTrack(Drp1Image,SelectedTracks,Drp1Rad,PixelSize,StartFrame,EndFrame, CircleMask)
% measure Drp1 intensity only on Drp1 punctae within the circle marking the
% constriction site
%
% Inputs:
% Drp1Image = path to Drp1 channel tif stack
% SelectedTracks = cell with different tracks, their ID, location and frame
% Drp1Rad = radius of Drp1 circle in nanometers 
% PixelSize = in [nm]
% StartFrame = 
% EndFrame
% 
% Outputs:
% TracksinFrames = matrix indicating which tracks are found in which frames
% Drp1Mask = binary image showing which regions are considered for
%            measuring the Drp1 intensity
% ROIsumfinal = measured intensity in Drp1Mask region for each frame


R0=Drp1Rad;
theta=0:.01:2*pi;

TracksinFrames=zeros(EndFrame,1);

% DonutMask=cell(1,size(TracksinFrames,1));

% %     for t=1:size(SelectedTracks,1)
% %         for f=SelectedTracks{t,1}(1,4):SelectedTracks{t,1}(end,4)
% %             if f>=StartFrame && f<=EndFrame
% %                 x0D(t,f)=SelectedTracks{t,1}(find(SelectedTracks{t,1}(:,4)==f),2)/PixelSize;
% %                 y0D(t,f)=SelectedTracks{t,1}(find(SelectedTracks{t,1}(:,4)==f),3)/PixelSize; 
% % 
% %                 xiDo{t,f}=2*R0*cos(theta)/PixelSize+x0D(t,f);
% %                 yiDo{t,f}=2*R0*sin(theta)/PixelSize+y0D(t,f);
% % 
% %                 imshow(imread(Drp1Image,f), [])
% %                 DrpCircleDo{t,f}=line(xiDo{t,f},yiDo{t,f}, 'LineWidth', 1.5, 'Color', [0 0.7 0] );
% %                 drawnow;
% %                 %mask making
% %                 roimaskDo{t,f}=poly2mask(xiDo{t,f},yiDo{t,f}, size(imread(Drp1Image, f),1), size(imread(Drp1Image, f),2));
% %                 
% %                 xi{t,f}=R0*cos(theta)/PixelSize+x0D(t,f);
% %                 yi{t,f}=R0*sin(theta)/PixelSize+y0D(t,f);
% % 
% %                 imshow(imread(Drp1Image,f), [])
% %                 DrpCircleDi{t,f}=line(xi{t,f},yi{t,f}, 'LineWidth', 1.5, 'Color', [0.7 0 0] );
% %                 drawnow;
% %                 %mask making
% %                 roimaskDi{t,f}=poly2mask(xi{t,f},yi{t,f}, size(imread(Drp1Image, f),1), size(imread(Drp1Image, f),2));
% %                 roimaskD{t,f}=roimaskDo{t,f}-roimaskDi{t,f};
% %                 
% %                 InMaskD{t,f}=find(roimaskD{t,f});
% %                 CurrentImageD{t,f}=imread(Drp1Image,f);
% %                 CurrentImageD{t,f}(CurrentImageD{t,f}<0)=0;
% %                 
% %                 ROImeanD(t,f)=mean(CurrentImageD{t,f}(InMaskD{t,f}));
% %             endDrp1Signal_BleachCorrected
% %     end
% % end
x0=zeros(size(SelectedTracks,1),EndFrame);
y0=zeros(size(SelectedTracks,1),EndFrame);
xi=cell(size(SelectedTracks,1),EndFrame);
yi=cell(size(SelectedTracks,1),EndFrame);
DrpCircle=cell(size(SelectedTracks,1),EndFrame);
roimask=cell(size(SelectedTracks,1),EndFrame);
InMask=cell(size(SelectedTracks,1),EndFrame);
CurrentImage=cell(size(SelectedTracks,1),EndFrame);
    for t=1:size(SelectedTracks,1)
        for f=SelectedTracks{t,1}(1,4):SelectedTracks{t,1}(end,4)
            if f>=StartFrame && f<=EndFrame
                
                x0(t,f)=SelectedTracks{t,1}(find(SelectedTracks{t,1}(:,4)==f),2)/PixelSize;
                y0(t,f)=SelectedTracks{t,1}(find(SelectedTracks{t,1}(:,4)==f),3)/PixelSize; 
                
%                 xi{t,f}=R0*cos(theta)/PixelSize+x0D(t,f);
%                 yi{t,f}=R0*sin(theta)/PixelSize+y0D(t,f);
%                 


                xi{t,f}=R0*cos(theta)/PixelSize+x0(t,f);
                yi{t,f}=R0*sin(theta)/PixelSize+y0(t,f);
                
                imshow(imread(Drp1Image,f), [])
                DrpCircle{t,f}=line(xi{t,f},yi{t,f}, 'LineWidth', 1.5, 'Color', [0.7 0 0] );
                drawnow;
                %mask making
                roimask{t,f}=poly2mask(xi{t,f},yi{t,f}, size(imread(Drp1Image, f),1), size(imread(Drp1Image, f),2));
    
                %apply mask to image
                InMask{t,f}=find(roimask{t,f});
                CurrentImage{t,f}=imread(Drp1Image,f);
                minCI=min(min(CurrentImage{t,f}));
%                 CurrentImage{t,f}=(CurrentImage{t,f})-minCI;
                
%                 ROImean(t,f)=mean(CurrentImage{t,f}(InMask{t,f}));
%                 ROIsum(t,f)=sum(sum(CurrentImage{t,f}(InMask{t,f})-ROImeanD(t,f)));
                
                TracksinFrames(f,length(nonzeros(TracksinFrames(f,:)))+1)=t;

            end
        end
    end
    
    Background=cell(1,size(TracksinFrames,1));
    AB=zeros(1,size(TracksinFrames,1));
    AT=zeros(1,size(TracksinFrames,1));
    InMaskT=cell(1,size(TracksinFrames,1));
    InMaskB=cell(1,size(TracksinFrames,1));
    IB=zeros(1,size(TracksinFrames,1));
    IT=zeros(1,size(TracksinFrames,1));
    I=zeros(1,size(TracksinFrames,1));
    for f=1:size(TracksinFrames,1)
        Background{f}=CircleMask{f};

        if (length(nonzeros(TracksinFrames(f,:)))>0)
            for i=1:length(nonzeros(TracksinFrames(f,:)))
                Background{f}=Background{f}-(roimask{TracksinFrames(f,i),f});
            end
        end
      
        Background{f}(Background{f}<=0)=0;
        Background{f}(Background{f}>0)=1;
        AB(f)=sum(sum(Background{f}));
        AT(f)=sum(sum(CircleMask{f}));
        
        InMaskT{1,f}=find(CircleMask{f});
        InMaskB{1,f}=find(Background{f});
        CurrentImage{1,f}=imread(Drp1Image,f);
        minCI=min(min(CurrentImage{1,f}));
%         CurrentImage{1,f}=(CurrentImage{1,f})-minCI;
        
%         ROImean(1,f)=mean(CurrentImage{1,f}(InMask{1,f})-ROImeanD(t,f));
%         ROIsumfinal(1,f)=sum(sum(CurrentImage{1,f}(InMask{1,f})-ROImeanD(t,f)));
        
        IB(f)=sum(sum(CurrentImage{1,f}(InMaskB{1,f})))*(AT(f)/AB(f));
        IT(f)=sum(sum(CurrentImage{1,f}(InMaskT{1,f})));
        I(f)=IT(f)-IB(f);
    end