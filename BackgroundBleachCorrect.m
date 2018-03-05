function [ Drp1Signal_BleachBkgrndCorrected, BleachSlope,BleachYInter, SubDrp1, AreaBox] = BackgroundBleachCorrect(StartFrame, EndFrame, DrpImage, MitoImage, measurement, frameStamp,  CircleMask, R0c, InMask)
%bleach correct the Drp1 signal after applying a background subtraction
ImD=zeros([size(imread(DrpImage,StartFrame)) (EndFrame)]);
ImM=zeros([size(imread(DrpImage,StartFrame)) (EndFrame)]);
for f=StartFrame:EndFrame
   IMD=imread(DrpImage,f);
   IMD(IMD<0)=0;
   ImD(:,:,f)=IMD;
   IMM=imread(MitoImage,f);
   IMM(IMM<0)=0;
   ImM(:,:,f)=IMM;
end
clear IMD; clear IMM;
MIP_Drp1=max(ImD, [],3);
MIP_mito=max(ImM, [],3);
DualColorMIP=imfuse(MIP_mito,(MIP_Drp1));

imshow(DualColorMIP,[])
hold on;
title('Draw a rectangle around background')

%user defined background

% [Xedge, Yedge]=getpts(1);
% BackgroundCircle=imellipse(gca, [Xedge, Yedge, 2*R0c,2*R0c]);
BackgroundPoly=impoly;
BackgroundPos=BackgroundPoly.getPosition;
 BackgroundPos1(:,1)=[BackgroundPos(:,1); BackgroundPos(1,1)];
 BackgroundPos1(:,2)=[BackgroundPos(:,2); BackgroundPos(1,2)];

close;
 mask(:,:)=poly2mask(BackgroundPos1(:,1),BackgroundPos1(:,2), size(ImD(:,:,1), 1), size(ImD(:,:,1),2));
% mask(:,:)=poly2mask(BackgroundPos1(:,1),BackgroundPos1(:,2), 173,173);

SubDrp1=zeros(size(ImD));
RawBackground=zeros(1,(EndFrame));
RawBackgroundAv=zeros(1,(EndFrame));
for t=StartFrame:EndFrame
% mask(:,:,t)=poly2mask(BackgroundPos1(:,1),BackgroundPos1(:,2), size(ImD(:,:,t), 1), size(ImD(:,:,t),2));
 SubDrp1(:,:,t)=ImD(:,:,t).*mask(:,:);
%   SDrp(SDrp<0)=0;
%   SubDrp1(:,:,t)=SDrp(:,:,t);
%   SubDrp1(:,:,t)=imcrop(ImD(:,:,t),BackgroundPos);
  RawBackground(t)=sum(sum(SubDrp1(:,:,t)));
  RawBackgroundAv(t)=mean(mean(SubDrp1(:,:,t)));
end

%AreaCMask=pi*(R0c)^2; %area of CircleMask around which signal is measured
AreaBox=polyarea(BackgroundPos1(:,1), BackgroundPos1(:,2));%size(SubDrp1(:,:,1),1)*size(SubDrp1(:,:,1),2);%BackgroundPos(3)*BackgroundPos(4); %area of box over which background was estimated
AreaCMask=size(InMask{1,1},1);

Bkgrnd=zeros(1,(EndFrame));
DrpConsite_BkgrnCorrected=zeros(1,(EndFrame));
for u=StartFrame:EndFrame
   Bkgrnd(u)=(AreaCMask/AreaBox)*RawBackground(u);
    DrpConsite_BkgrnCorrected(u)=measurement(u)- Bkgrnd(u);
end



plot(StartFrame:EndFrame, RawBackground, '-')
hold on
plot(StartFrame:EndFrame, RawBackground, 'g+')
title('Background signal versus time');
close;

LinearBleachFit=polyfit(StartFrame:EndFrame, RawBackground,1);
BleachSlope=LinearBleachFit(1,1);
BleachYInter=LinearBleachFit(1,2);

%apply to Drp1 signal measurement
Eff=zeros(1,(EndFrame));
Drp1Signal_BleachBkgrndCorrected=zeros(1,(EndFrame));
for i=frameStamp(1):frameStamp(end)
    Eff(i)=((BleachSlope*i)+BleachYInter)./BleachYInter; %fluorophore efficiency
    Drp1Signal_BleachBkgrndCorrected(i)=DrpConsite_BkgrnCorrected(i)./Eff(i);
end

plot(measurement, 'g+')
hold on;
plot(measurement, 'b-')
plot(Bkgrnd, 'm+');
plot(Bkgrnd, 'r-');
title('Raw data (Green) and background (Pink)')

figure;
plot(Drp1Signal_BleachBkgrndCorrected, 'b*')
hold on
plot(Drp1Signal_BleachBkgrndCorrected, 'k-')
title('background and bleach corrected Drp1 signal')


    
end


