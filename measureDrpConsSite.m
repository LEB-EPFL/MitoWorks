function [ CircleMask,Drp1IntensityConsSite,R0c, InMask] = measureDrpConsSite(ContnextXsnake, ContnextYsnake, IDmin, subNormals,TimeStamp, StartFrame, EndFrame, Drp1Image, Pixel)
%this function draws a circle at the min diameter of each frame and
%measures the Drp1 flux within this circle
%Drp1IntensityConsSite is a matrix with as many entries as frames analyzed

disp('Drp1 measurement radius');
Drp1RadCon=askDrp1MeasureRadius;

R0c=Drp1RadCon./Pixel;
q=0:.01:2*pi;
X_MinConstriction=zeros((EndFrame),1);
Y_MinConstriction=zeros((EndFrame),1);
Xcircle=cell(1,(EndFrame));
Ycircle=cell(1,(EndFrame));
Circle=cell(1,(EndFrame));
CircleMask=cell(1,(EndFrame));
InMask=cell(1,(EndFrame));
CurrentImage=cell(1,(EndFrame));
ROImean=cell(1,(EndFrame));
ROIsum=zeros(1,(EndFrame));
%find centre location of min diameter (at centre of backbone/centre line) 
for f=StartFrame:EndFrame
    
    X_MinConstriction(f,1)=mean(subNormals{1,f}(IDmin(f),1:2:3));
    Y_MinConstriction(f,1)=mean(subNormals{1,f}(IDmin(f),2:2:4));
    
    Xcircle{1,f}=R0c*cos(q)+X_MinConstriction(f,1);
    Ycircle{1,f}=R0c*sin(q)+Y_MinConstriction(f,1);
    
    figure;
    imshow(imread(Drp1Image, f), [])
    hold on
    plot(ContnextXsnake{1,f},ContnextYsnake{1,f}, 'LineWidth', 0.5, 'Color', [.7 0 0.4 ])
    Circle{1,f}=line(Xcircle{1,f}, Ycircle{1,f}, 'LineWidth', 1.5, 'Color', [0 0.5 .9] );
    
     %mask making
    CircleMask{1,f}=poly2mask(Xcircle{1,f}, Ycircle{1,f}, size(imread(Drp1Image,f),1), size(imread(Drp1Image, f),2));
    
    %apply mask to image
    InMask{1,f}=find(CircleMask{1,f});
    CurrentImage{1,f}=imread(Drp1Image,f);
    
    CurrentImage{1,f}(CurrentImage{1,f}<0)=0;
    
    %measure stats
    
    ROImean{1,f}=mean(mean(CurrentImage{1,f}(InMask{1,f})));
    ROIsum(1,f)=sum(sum(CurrentImage{1,f}(InMask{1,f})));
    
end



% ROIsum=cell2mat(ROIsum);
ROIsum(~any(ROIsum,2), :)=[];
Drp1IntensityConsSite=ROIsum;

close all;
time=TimeStamp;
time(~any(time(1:end-1),2), :)=[]; %delete rows with only zeros
figure
title ('Drp1 flux at constriction site versus time')

plot(StartFrame:EndFrame,ROIsum(1,StartFrame:EndFrame), '*' )
xlabel('Time before fission (s)')
ylabel ('Drp1 flux')

