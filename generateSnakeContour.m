%script that generates a mitochondrial contour, calculates
%envelope curvature along contour and measures distnaces across mito
%according to the intensity histograms

function [ContnextXsnake, ContnextYsnake, pos ]=generateSnakeContour(SegmentedOriginalImage,CropData,StartFrame,EndFrame,px)

% fname = '151128_comp119_event1GREEN-1.tif' % name of mitochondria probability map
% CropData = '119fullMovie-1.tif' %name of mito raw data-- this MUST be the same size (dimensions, no of frames) as the probability map
% nframe = 2% number of frames in movie before break
%StartFrame=1; %frame to start analysis
% px = 30 %pixel size in nm
% S = 28%smoothing S*10 (in nm)
nVisIn = EndFrame; %number of frames to be visually checked 
% meshName= '119MitoMesh.mat'

% addpath('\\ipsbsrv4.epfl.ch\IPSB\LEB\Shared\LinaTatjana\contours_Scripts\MicrobeTracker 0.936'); %change to location of MicrobeTracker

%import the raw data tiff stack (16 bit)

%import the prob map tiff stack (16 bit)

%crop ROI 
h = waitbar(0,'Generating snake contour...');
I=zeros([size(SegmentedOriginalImage{StartFrame}) (EndFrame)]);
Ibl=zeros([size(SegmentedOriginalImage{StartFrame}) (EndFrame)]);
Image=zeros(3*[size(SegmentedOriginalImage{StartFrame}) (EndFrame)]);
I_filter=zeros(3*[size(SegmentedOriginalImage{StartFrame}) (EndFrame)]);
for f=StartFrame:EndFrame
            waitbar((f-StartFrame)/(EndFrame-StartFrame));

   I(:,:,f)=SegmentedOriginalImage{f};
   Ibl(:,:,f)=imgaussfilt(I(:,:,f), 0.2);
   Image(:,:,f)=imresize(Ibl(:,:,f), 3);
   %imshow(image);
  % hold on
   %h=imfreehand();
   %binaryImage=h.createMask();
   %close
   %Inv(:,:,i)=imcomplement(binaryImage);
   %image(Inv(:,:,i))=255;
   I_filter(:,:,f)=Image(:,:,f);

end
close(h);
% % %ask user to specify experiment
% % choice = questdlg('Please select experiment', ...
% %     'Experiment', ...
% %     'Tom20','MitoTracker', 'MitoTracker');
% % % Handle response
% % switch choice
% %     case 'Tom20'
% %         exp=1;
% %         for f=StartFrame:EndFrame
% %             I_filter(:,:,f)=imfill(I_filter(:,:,f));
% %         end
% %     case 'MitoTracker'
% %         exp=2;
% % end

%active contouring using 'Chan-Vese' snake

% compose total image
ComposedImage=(zeros(size(I_filter(:,:,StartFrame))));
for f=StartFrame:EndFrame
    ComposedImage=max(ComposedImage,I_filter(:,:,f));
end
ComposedImageBW=imdilate(im2bw(ComposedImage),strel('ball',round(300/px),0));
% size(ComposedImage)
% imshow(ComposedImageBW);
% title('drag a box around mito')
% box=impoly;
% pos=box.getPosition();
pos=[0 0 0 0];
% for f=StartFrame:EndFrame
%    Icrop(:,:,f)=imcrop(I_filter(:,:,f),pos); 
% end

% boundary=box.createMask();
boundary=ComposedImageBW;
% mask=uint16(ones(size(Icrop(:,:,StartFrame))));
% mask=+boundary;

% box=imrect;
% boundary=createMask(box);
% mask=uint16(zeros(size(I_filter{StartFrame})));
% mask=+boundary;
% size(mask)
% size(I_filter{1})
S=askSmoothingFactor();
%         H= waitbar(0,'Generating initial contour');
parfor f=StartFrame:EndFrame
    fprintf('|');
   bw(:,:,f)=activecontour(I_filter(:,:,f),boundary,1000, 'Chan-Vese','SmoothFactor', 1.5);

   %generate approximate contour
  conture = imcontour(bw(:,:,f), 1)
   ApC{:,f}=conture(:,2:end);
%    close;
% ApC{:,f}=imcontour(bw(:,:,f), 1);
%    ApC{:,f}=ApC{1,f}(:,2:end);
   
   %smooth the contour
   smContour{:,f}=lssmooth(ApC{:,f}',0.5*str2num(S{1}));
   
   %resize the smooth contour
   Contour_f{:,f}=[smContour{:,f}./3; smContour{:,f}(1,:)./3];
%    Contour_f{:,f}(end+1,:)=Contour_f{:,f}(1,:)
   
   %smooth again
   Contour_final{:,f}=lssmooth(Contour_f{:,f}, str2num(S{1}));
  
end
% close(H);

choice1=questdlg('Are you ready to inspect the contour?',...
    'Contour check',...
    'Go!', 'Go!' );
%Handle button click and visually inspect contour over SIM data
switch choice1
    case 'Go!'
     exp=1;
     
ContnextXsnake=cell(1,(nVisIn-StartFrame+1));
ContnextYsnake=cell(1,(nVisIn-StartFrame+1));
for ii=StartFrame:nVisIn 

  figure;
%    axis equal
   imshow(imread(CropData,ii), []);
   hold on;
   plot(Contour_final{1,ii}(:,1),Contour_final{1,ii}(:,2), 'y-', 'linewidth', 1.8);
   drawnow;
   str = sprintf('frame %d ', ii);
   title(str);
   ContnextXsnake{ii}=(Contour_final{1,ii}(:,1)); %+pos(1)/3);%.*px;
   ContnextYsnake{ii}=(Contour_final{1,ii}(:,2));%+pos(2)/3);%.*px;
%    ContnextXsnake{ii}(end+1)=ContnextXsnake{ii}(1);
%    ContnextYsnake{ii}(end+1)=ContnextYsnake{ii}(1);
   
end
h=msgbox('Press any key to continue once you have checked all frames and noted mistakes (if any)')
pause
% close h;
end
close all; 

choice2=questdlg('Does the contour look ok for every frame? Warning: a contour that does not converge will cause problems if you proceed',...
    'Contour check on each frame',...
    'I want to re-run one or more frames', 'Proceed, everything looks OK', 'Proceed, everything looks OK');
switch choice2
  case 'I want to re-run one or more frames'
        exp=1;
       FramesToCorrectSnake=askFrames();    
       XsnakeCorrect=cell(1,(EndFrame));
       YsnakeCorrect=cell(1,(EndFrame));
     for fr=1:length(FramesToCorrectSnake);
           [XsnakeCorrect{FramesToCorrectSnake(fr)}, YsnakeCorrect{FramesToCorrectSnake(fr)}]=CorrectSnake(SegmentedOriginalImage{FramesToCorrectSnake(fr)},imread(CropData, FramesToCorrectSnake(fr)), FramesToCorrectSnake(fr), px, S);
   
       ContnextXsnake{1,FramesToCorrectSnake(fr)}=XsnakeCorrect{FramesToCorrectSnake(fr)}{1,FramesToCorrectSnake(fr)};
       ContnextYsnake{1,FramesToCorrectSnake(fr)}=YsnakeCorrect{FramesToCorrectSnake(fr)}{1,FramesToCorrectSnake(fr)};
       %plot it
       figure;
       imshow(imread(CropData,FramesToCorrectSnake(fr)), []);
       hold on
 plot(ContnextXsnake{1,FramesToCorrectSnake(fr)},ContnextYsnake{1,FramesToCorrectSnake(fr)}, 'y-');
       
     end

       
    case 'Proceed, everything looks OK'
        
end

end

   