function [XsnakeCorrect,YsnakeCorrect ]=CorrectSnake(SegmentedOriginalImage,CropData,FramesToCorrectSnake,px, S )

% fname = '151128_comp119_event1GREEN-1.tif' % name of mitochondria probability map
% CropData = '119fullMovie-1.tif' %name of mito raw data-- this MUST be the same size (dimensions, no of frames) as the probability map
% nframe = 2% number of frames in movie before break
%StartFrame=1; %frame to start analysis
% px = 30 %pixel size in nm
% S = 28%smoothing S*10 (in nm)
% nVisIn = EndFrame; %number of frames to be visually checked 
% meshName= '119MitoMesh.mat'

% addpath('\\ipsbsrv4.epfl.ch\IPSB\LEB\Shared\LinaTatjana\contours_Scripts\MicrobeTracker 0.936'); %change to location of MicrobeTracker

%import the raw data tiff stack (16 bit)

%import the prob map tiff stack (16 bit)

%crop ROI 
h = waitbar(0,'Generating snake contour...');

g=FramesToCorrectSnake;

%    I(:,:,g)=SegmentedOriginalImage{g};
      I(:,:,g)=SegmentedOriginalImage;

   Ibl(:,:,g)=imgaussfilt(I(:,:,g), 0.2);
   Image(:,:,g)=imresize(Ibl(:,:,g), 3);
   %imshow(image);
  % hold on
   %h=imfreehand();
   %binaryImage=h.createMask();
   %close
   %Inv(:,:,i)=imcomplement(binaryImage);
   %image(Inv(:,:,i))=255;
   I_filter(:,:,g)=Image(:,:,g);


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
ComposedImage=(zeros(size(I_filter(:,:,g))));

    ComposedImage=max(ComposedImage,I_filter(:,:,g));

ComposedImageBW=imdilate(im2bw(ComposedImage),strel('ball',round(300/px),0));
% ComposedImageBW1=imfill(ComposedImageBW);
% size(ComposedImage)
% imshow(ComposedImageBW1, []);
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

%S=askSmoothingFactor();
%         H= waitbar(0,'Generating initial contour');


% for f=StartFrame:EndFrame
%     f=FramesToCorrectSnake(f);

%   disp(sprintf('Frame %d / %d',f-StartFrame+1,EndFrame-StartFrame+1));
   bw(:,:,g)=activecontour(I_filter(:,:,g),boundary,1000, 'Chan-Vese','SmoothFactor', 1.5);

   %generate approximate contour
  conture = imcontour(bw(:,:,g), 1);
   ApC{:,g}=conture(:,2:end);
 close;
 ApC{:,g}(end+1,:)=ApC{:,g}(1,:);
% ApC{:,f}=imcontour(bw(:,:,f), 1);
%    ApC{:,f}=ApC{1,f}(:,2:end);
   
   %smooth the contour
   smContour{:,g}=lssmooth(ApC{:,g}',0.5*str2num(S{1}));
   
   %resize the smooth contour
   Contour_f{:,g}=smContour{:,g}./3;
   
   %smooth again
   Contour_final{:,g}=lssmooth(Contour_f{:,g}, str2num(S{1}));
  XsnakeCorrect{g}= (Contour_final{1,g}(:,1));
   
   YsnakeCorrect{g}=(Contour_final{1,g}(:,2));

end






