function [ Drp1IntensityFinalTrack] = measureDrpSignal_finalTrack(Drp1Image, DrpTrackCons, Drp1Rad, PixelSize, TimeStamp,  TimeRes,StartFrame, EndFrame )

%Drp1IntIntensityFinalTrack is a matrix, where each entry is the Drp1 flux
%(total intensity/area) within a user defined circle centered at
%thecentroid of each Drp1 spot; there are as many columns as there are
%frames


%draw a circle around each Drp1 spot, measure intensity inside this ROI

R0=Drp1Rad;
theta=0:.01:2*pi;


x0D=zeros(1,EndFrame);
y0D=zeros(1,EndFrame);
xiDo=cell(1,EndFrame);
yiDo=cell(1,EndFrame);
DrpCircleDo=cell(1,EndFrame);
roimaskDo=cell(1,EndFrame);
xi=cell(1,EndFrame);
yi=cell(1,EndFrame);
DrpCircleDi=cell(1,EndFrame);
roimaskDi=cell(1,EndFrame);
roimaskD=cell(1,EndFrame);
InMaskD=cell(1,EndFrame);
CurrentImageD=cell(1,EndFrame);
ROImeanD=zeros(1,EndFrame);
Circle=cell(1,EndFrame);
roimask=cell(1,EndFrame);
InMask=cell(1,EndFrame);
CurrentImage=cell(1,EndFrame);
ROImean=zeros(1,EndFrame);
ROIsum=zeros(1,EndFrame);
for f=DrpTrackCons(1,4):DrpTrackCons(end,4)
    if f>=StartFrame && f<=EndFrame
        x0D(1,f)=DrpTrackCons(find(DrpTrackCons(:,4)==f),2)./PixelSize;
        y0D(1,f)=DrpTrackCons(find(DrpTrackCons(:,4)==f),3)./PixelSize; 

        xiDo{1,f}=2*R0*cos(theta)/PixelSize+x0D(1,f);
        yiDo{1,f}=2*R0*sin(theta)/PixelSize+y0D(1,f);

        imshow(imread(Drp1Image,f), [])
        DrpCircleDo{1,f}=line(xiDo{1,f},yiDo{1,f}, 'LineWidth', 1.5, 'Color', [0 0.7 0] );
        drawnow;
        hold on
        %mask making
        roimaskDo{1,f}=poly2mask(xiDo{1,f},yiDo{1,f}, size(imread(Drp1Image, f),1), size(imread(Drp1Image, f),2));

        xi{1,f}=R0*cos(theta)/PixelSize+x0D(1,f);
        yi{1,f}=R0*sin(theta)/PixelSize+y0D(1,f);

%         imshow(imread(Drp1Image,f), [])
        DrpCircleDi{1,f}=line(xi{1,f},yi{1,f}, 'LineWidth', 1.5, 'Color', [0 0.7 0] );
        drawnow;
        hold off
        %mask making
        roimaskDi{1,f}=poly2mask(xi{1,f},yi{1,f}, size(imread(Drp1Image, f),1), size(imread(Drp1Image, f),2));
        roimaskD{1,f}=roimaskDo{1,f}-roimaskDi{1,f};

        InMaskD{1,f}=find(roimaskD{1,f});
        CurrentImageD{1,f}=imread(Drp1Image,f);
        minCID=min(min(CurrentImageD{1,f}));
%         CurrentImageD{1,f}=(CurrentImageD{1,f})-minCID;

        ROImeanD(f)=mean(mean(CurrentImageD{1,f}(InMaskD{1,f})));
    end
end

%make a circular ROI at each frame, define a mask with it and measure
%intensity within this mask
for k=DrpTrackCons(1,4):DrpTrackCons(end,4)
        if k>=StartFrame && k<=EndFrame

    xi{1,k}=(R0*cos(theta)/PixelSize+x0D(k));
    yi{1,k}=(R0*sin(theta)/PixelSize+y0D(k));
%     figure;
    imshow(imread(Drp1Image,k), [])
    hold on
    Circle{1,k}=line(xi{1,k}, yi{1,k}, 'LineWidth', 1.5, 'Color', [0.7 0 0] );
    drawnow
    hold off
    %mask making
    roimask{1,k}=poly2mask(xi{1,k}, yi{1,k}, size(imread(Drp1Image, k),1), size(imread(Drp1Image, k),2));
    
    %apply mask to image
    InMask{1,k}=find(roimask{1,k});
    CurrentImage{1,k}=imread(Drp1Image,k);
    minCI=min(min(CurrentImage{1,k}));
%     CurrentImage{1,k}=(CurrentImage{1,k})-minCI;
    
    %measure stats
    
    ROImean(1,k)=mean(mean(CurrentImage{1,k}(InMask{1,k})-ROImeanD(k)));
    ROIsum(1,k)=sum(sum(CurrentImage{1,k}(InMask{1,k})-ROImeanD(k)));
    
% ROIstd{1,k}=std(imread(Drp1Image,DrpTrackCons(k,4))(InMask{1,k}));

        end
end

Drp1IntensityFinalTrack=ROIsum;





close all;
time=zeros(1,length(size(Drp1IntensityFinalTrack, 2)));
timeSHIFT=zeros(1,length(size(Drp1IntensityFinalTrack, 2)));
for i=1:size(Drp1IntensityFinalTrack, 2);
time(i)=i*TimeRes;
timeSHIFT(i)=time(i)-size(Drp1IntensityFinalTrack, 2)*TimeRes;
end
%time(~any(time(1:EndFrame-1,:),2), :)=[]; %delete rows with only zeros
figure
title ('Drp1 flux versus time')

plot(timeSHIFT',Drp1IntensityFinalTrack, '*' )
xlabel('Time before fission (s)')
ylabel ('Drp1 flux')

%%TO DO: plot bleach corrected versus non-bleach corrected signal