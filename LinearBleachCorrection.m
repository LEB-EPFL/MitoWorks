function [ Drp1Signal_BleachCorrected] = LinearBleachCorrection(measurement, frameStamp, BleachSlope,BleachYInter)
%bleach correct the Drp1 signal

% for f=StartFrame:EndFrame
%    ImD(:,:,f)=imread(DrpImage,f);
%    ImM(:,:,f)=imread(MitoImage,f);
% end
% MIP_Drp1=max(ImD, [],3);
% MIP_mito=max(ImM, [],3);
% DualColorMIP=imfuse(MIP_mito,imadjust(MIP_Drp1));
% 
% imshow(DualColorMIP,[])
% hold on;
% title('draw box where you want to measure background bleaching')
% 
% %user defined background
% BackgroundBox=imrect;
% BackgroundPos=BackgroundBox.getPosition();
% close;
% 
% for t=StartFrame:EndFrame
%   SubDrp1(:,:,t)=imcrop(ImD(:,:,t),BackgroundPos);
%   RawBackground(t)=sum(sum(SubDrp1(:,:,t)));
%   RawBackgroundAv(t)=mean(mean(SubDrp1(:,:,t)));
% end
% plot(StartFrame:EndFrame, RawBackground, '-')
% hold on
% 
% plot(StartFrame:EndFrame, RawBackground, 'g+')
% 
% LinearBleachFit=polyfit(StartFrame:EndFrame, RawBackground,1);
% Slope=LinearBleachFit(1,1);
% YInter=LinearBleachFit(1,2);

%apply to Drp1 signal measurement
Drp1Signal_BleachCorrected=zeros(1,frameStamp(end));
Eff=zeros(1,frameStamp(end));
for i=frameStamp(1):frameStamp(end)
    Eff(i)=((BleachSlope*i)+BleachYInter)./BleachYInter; %fluorophore efficiency
%     Drp1Signal_BleachCorrected(i)= measurement(i+1-(frameStamp(1)))./Eff(i);
Drp1Signal_BleachCorrected(i)= measurement(i)./Eff(i);
end

plot(Drp1Signal_BleachCorrected(frameStamp(1):frameStamp(end)), 'b*')
hold on
plot(Drp1Signal_BleachCorrected(frameStamp(1):frameStamp(end)), 'k-')
title('bleach corrected Drp1 signal')


    
end


