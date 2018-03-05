function [ ContourSnake ] = generateFinalContours(SegmentedMitoIDImage, OriginalImage, StartFrame, EndFrame)
%generateFinalContours produces a final mitochondria contour and 2 contours for
%comparison based on the SIM image

%   Detailed explanation goes here
%script that generates a mitochondrial contour, calculates
%envelope curvature along contour and measures distances across mito
%according to the intensity histograms

% % % clear all
% % % close all
% % % addpath('F:\Data-E\Olympus\2015\CurvatureScripts\ActiveContour\TIFFStack\private')
% % % fname = 'Probability maps-newNOblur.tif'; % name of mitochondria probability map
% % % CropData = '119OriginalMovie1-23.tif' %name of mito raw data-- this MUST be the same size (dimensions, no of frames) as the probability map
% % % nframe = 21;% number of frames in movie before break
% % % start_frame=17; %frame to start_frame analysis
% % % px = 30 %pixel size in nm
% % % S = 28%smoothing S*10 (in nm)
% % % nVisIn = nframe %number of frames to be visually checked 
% % % meshName= '119newOrderMitoMesh.mat'
% % % meshSpacing=1;


% addpath('\\ipsbsrv4.epfl.ch\IPSB\LEB\Shared\LinaTatjana\contours_Scripts\MicrobeTracker 0.936'); %change to location of MicrobeTracker

%import the raw data tiff stack (16 bit)
% % cData=TIFFStack(CropData);
% % 
% % %import the prob map tiff stack (16 bit)
% % Stack=TIFFStack(fname);
% % 
% % %crop ROI 


for i=StartFrame:EndFrame
   SegmentedMitoIDImage{f}=im2uint16(SegmentedMitoIDImage{f});
   I(:,:,f)=SegmentedMitoIDImage{f};
   Ibl(:,:,f)=imgaussfilt(I(:,:,f), 1);
   Image(:,:,f)=imresize(Ibl(:,:,f), 3);
   I_filter(:,:,f)=Image(:,:,f);
    
%    I(:,:,i)=Stack(:,:,i);
%    Ibl(:,:,i)=imgaussfilt(I(:,:,i), 1);
%    image=imresize(Ibl(:,:,i), 3);
%    imshow(image);
%    hold on
%    h=imfreehand();
%    binaryImage=h.createMask();
%    close
%    Inv(:,:,i)=imcomplement(binaryImage);
%    image(Inv(:,:,i))=255;
%    I_filter(:,:,i)=image;

end

%ask user to specify experiment
% % choice = questdlg('Please select experiment', ...
% %     'Experiment', ...
% %     'Tom20','MitoTracker', 'MitoTracker');
% % % Handle response
% % switch choice
% %     case 'Tom20'
% %         exp=1;
% %         for k=1:nframe
% %             I_filter(:,:,k)=imfill(I_filter(:,:,k));
% %         end
% %     case 'MitoTracker'
% %         exp=2;
% % end

%active contouring using 'Chan-Vese' snake
% % mip=max(I_filter, [], 3);
% % imshow(mip);
% % box=imrect;
% % boundary=box.createMask();
% % mask=+boundary;

ComposedImage=uint16(zeros(size(I_filter(:,:,StartFrame))));
for f=StartFrame:EndFrame
    ComposedImage=max(ComposedImage,I_filter(:,:,f));
end

% size(ComposedImage)
imshow(ComposedImage);

box=imrect;
pos=box.getPosition();
for f=StartFrame:EndFrame
   Icrop(:,:,f)=imcrop(I_filter(:,:,f),pos); 
end

boundary=box.createMask();
mask=uint16(ones(size(Icrop(:,:,StartFrame))));

S=askSmoothingFactor();

for j=StartFrame:EndFrame
   bw{:,j}=activecontour(Icrop(:,:,j),mask, 1500, 'Chan-Vese', 'SmoothFactor', 1.5);
   
   %generate approximate contour
   ApC{:,j}=imcontour(bw{1,j}, 1);
   close;
   ApC{:,j}=ApC{1,j}(:,2:end);
   
   %smooth the contour
   smContour{:,j}=lssmooth(ApC{:,j}',S );
   
   %resize the smooth contour
   Contour_f{:,j}=smContour{:,j}./3;
   
   %smooth again
   ContourSnake{:,j}=lssmooth(Contour_f{:,j}, S)
end

% %visually inspect contour over SIM data
% for ii=start_frame:nVisIn 
%    h=figure;
%    axis equal
%    imshow(cData(:,:,ii), []);
%    hold on;
%    plot(Contour_final{1,ii}(:,1),Contour_final{1,ii}(:,2), 'y-', 'linewidth', 1.8);
% end

% %calculate the curvature of entire contour for every frame
% 
% for jj=start_frame:nframe
%     
%    ContnextX{:,jj}=Contour_final{:,jj}(:,1);
%    ContnextY{:,jj}=Contour_final{:,jj}(:,2);
%    v{jj}=[ContnextX{:,jj}, ContnextY{:,jj}];
%    Lines{jj}=[(1:(size(v{1,jj},1)-2))' (3:size(v{1,jj},1))'];% the curvature at a single point uses 3 neighbouring points (7 total: itself and 3 on each side)
%     
%   [k0{jj}]=LineCurvature2D(v{jj}, Lines{jj});
%   [k1{jj}]=[k0{jj}./px];
%     
% end

%ask user if normal lines across contour exist
choice = questdlg('Did you run MitoMesh.m before this session?', ...
    'Answer', ...
    'Yes','No', 'No');
% Handle response
switch choice
    case 'Yes'
        exp=1;
    MeshData=importdata(meshName);
    Norm=MeshData.subNormals;
    mesh=MeshData.mesh;
    
    for i=start_frame:nframe
    for j=1:length(Norm{i})
        
   Norm{:,i}(all(Norm{:,i}==0,2), :)=[];  
   
   TotDist{:,i}=sqrt((Norm{:,i}(:,1)-Norm{:,i}(:,3)).*(Norm{:,i}(:,1)-Norm{:,i}(:,3))+(Norm{:,i}(:,2)-Norm{:,i}(:,4)).*(Norm{:,i}(:,2)-Norm{:,i}(:,4)));
    end
    end
    
%     Norm=cellfun(@(x) x*px,Norm1,'un',0);

 for i=start_frame:nframe
    for j=1:length(Norm{i})
        
[cx{j,i}, cy{j,i}, c{j,i}, xi{j,i}, yi{j,i}]=improfile(cData(:,:,i), [Norm{1,i}(j,1), Norm{1,i}(j,3)], [Norm{1,i}(j,2), Norm{1,i}(j,4)], 'bilinear' );

%cx, cy are the spatial coordinates in pixels across the profile, c is the
%intensity of the profile, xi and yi are the inital points

%subtract background from intensity c{j,i} 
MinI{j,i}=min(c{j,i});
Int{j,i}=c{j,i}-MinI{j,i}; %background subtracted profile

%transform real space coordinates to distance along line and compute
%integrated intensity for each profile
Dx{j,i}=cx{j,i}-xi{j,i}(1,1);
Dy{j,i}=cy{j,i}-yi{j,i}(1,1);
dist{j,i}= sqrt((Dx{j,i}.*Dx{j,i})+(Dy{j,i}.*Dy{j,i})); %this is the x distance along x-axis
IntInt{1,i}(j,:)=trapz(dist{j,i}, Int{j,i}); %integrated intensity for each profile

%calculate truncated cone area 
radii{1,i}=sqrt((mesh{1,i}(:,1)-mesh{1,i}(:,3)).^2+(mesh{1,i}(:,2)-mesh{1,i}(:,4)).^2)/2;
for r=1:(length(mesh{1,i})-1)
subsegmentArea{1,i}(r,:)=(pi()*(radii{1,i}(r,1)+radii{1,i}(r+1,1))*sqrt((radii{1,i}(r,1)-radii{1,i}(r+1,1))^2+(meshSpacing)^2));
RatioR{1,i}(r, :)=radii{1,i}(r+1,1)/radii{1,i}(r,1);

% IntInt_av{1,i}(r,:)= (IntInt{1,i}(r,1)+IntInt{1,i}(r+1,1))/2;
% Int_Dis{1,i}(r,:)= IntInt_av{1,i}(r,:)./IntInt_av{1,i}(r,:);
end


% MitoIntInt_av(1,i)=sum(IntInt_av{1,i}); %sum integrated intensity for full mito
% %perform gaussian fit
% %fit_Gauss{j,i}=fit(dist{j,i}, Int{j,i}, 'gauss2');


% %find the FWHM points of intersection with x-axis
[width{j,i}, xLeft{j,i}, xRight{j,i},yLeft{j,i},yRight{j,i}]= getFWHM(dist{j,i},Int{j,i});
% [width{1,i}(j, :), xLeft{1,i}(j, :), xRight{1,i}(j,:),yLeft{1,i}(j, :),yRight{1,i}(j,:)]= getFWHM(dist{j,i},c{j,i});


% %transform back into real space
delta_x{j,i}=xi{j,i}(2)-xi{j,i}(1);
delta_y{j,i}=yi{j,i}(2)-yi{j,i}(1);
delta_yx{j,i}=(delta_y{j,i})./(delta_x{j,i});

delta{j,i}=[delta_x{j,i}, delta_y{j,i}];

if sign(delta{j,i}) == [1 1] %positve x and  y displacement
 
theta{j,i}= atan(delta_yx{j,i});
 
x_shift1{j,i}=xLeft{j,i}*abs(cos(theta{j,i}));
y_shift1{j,i}=xLeft{j,i}*abs(sin(theta{j,i}));

x_shift2{j,i}=xRight{j,i}*abs(cos(theta{j,i}));
y_shift2{j,i}=xRight{j,i}*abs(sin(theta{j,i}));
        
ConMito{1,i}(j, :) = [(xi{j,i}(1)+x_shift1{j,i}), (yi{j,i}(1)+y_shift1{j,i}), (xi{j,i}(1)+x_shift2{j,i}),(yi{j,i}(1)+y_shift2{j,i})];
%ConMito2{j,i} = [(xi{j,i}(1)+x_shift2{j,i}),(yi{j,i}(1)+y_shift2{j,i})];

end

 if sign(delta{j,i})== [-1 -1] %negative x and y displacement
 theta{j,i}= atan(delta_yx{j,i});
 alpha{j,i}=pi- theta{j,i};
 
x_shift1{j,i}=xLeft{j,i}*abs(cos(alpha{j,i}));
y_shift1{j,i}=xLeft{j,i}*abs(sin(alpha{j,i}));

x_shift2{j,i}=xRight{j,i}*abs(cos(alpha{j,i}));
y_shift2{j,i}=xRight{j,i}*abs(sin(alpha{j,i}));
      
ConMito{1,i}(j, :) = [(xi{j,i}(1)-x_shift1{j,i}), (yi{j,i}(1)-y_shift1{j,i}), (xi{j,i}(1)-x_shift2{j,i}),(yi{j,i}(1)-y_shift2{j,i})];
%ConMito2{j,i} = [(xi{j,i}(1)-x_shift2{j,i}),(yi{j,i}(1)-y_shift2{j,i})];
  end


 if sign(delta{j,i})== [1 -1] %positive x and negative y displacement
theta{j,i}= atan(delta_yx{j,i});
 
x_shift1{j,i}=xLeft{j,i}*abs(cos(theta{j,i}));
y_shift1{j,i}=xLeft{j,i}*abs(sin(theta{j,i}));

x_shift2{j,i}=xRight{j,i}*abs(cos(theta{j,i}));
y_shift2{j,i}=xRight{j,i}*abs(sin(theta{j,i}));
        
ConMito{1,i}(j, :) = [(xi{j,i}(1)+x_shift1{j,i}), (yi{j,i}(1)-y_shift1{j,i}), (xi{j,i}(1)+x_shift2{j,i}),(yi{j,i}(1)-y_shift2{j,i})];
%ConMito2{j,i} = [(xi{j,i}(1)+x_shift2{j,i}),(yi{j,i}(1)-y_shift2{j,i})];
 end

 if sign(delta{j,i})== [-1 1] %negative x and positive y displacement
        theta{j,i}= atan(delta_yx{j,i});
 
x_shift1{j,i}=xLeft{j,i}*abs(cos(theta{j,i}));
y_shift1{j,i}=xLeft{j,i}*abs(sin(theta{j,i}));

x_shift2{j,i}=xRight{j,i}*abs(cos(theta{j,i}));
y_shift2{j,i}=xRight{j,i}*abs(sin(theta{j,i}));

ConMito{1,i}(j, :) = [(xi{j,i}(1)-x_shift1{j,i}), (yi{j,i}(1)+y_shift1{j,i}), (xi{j,i}(1)-x_shift2{j,i}),(yi{j,i}(1)+y_shift2{j,i})];      
% ConMito{j,i} = [(xi{j,i}(1)-x_shift1{j,i}), (yi{j,i}(1)+y_shift1{j,i}), (xi{j,i}(1)-x_shift2{j,i}),(yi{j,i}(1)+y_shift2{j,i})];
% %ConMito2{j,i} = [(xi{j,i}(1)-x_shift2{j,i}),(yi{j,i}(1)+y_shift2{j,i})];
 end
 
 
    end
    %final smooth of FWHM contour
    
 %remove NaN values 
         
% [rN{1,i}, cN{1,i}]=find(isnan(ConMito{1,i}));
% sort(rN{1,i});
% rNshift{1,i}=circshift(rN{1,i}(:,1),1);
% row1{1,i}=abs(rNshift{1,i}(:,1)- rN{1,i}(:,1));
% rowNaN{1,i}=rN{1,i}(row1{1,i});

 %remove NaN values 
FindNaN= cellfun(@(x) (isnan(x)),ConMito, 'UniformOutput', false);
[rNaN{i}, cNaN{i}]=find(FindNaN{1,i}==1);
rNaNsort{1,i}=sort(rNaN{1,i});
rNshift{1,i}=circshift(rNaNsort{1,i},1);
[row1{1,i}, col1{1,i}]=find(abs(rNshift{1,i}(:,1)- rNaNsort{1,i}(:,1))==1);

OneNaN=any(FindNaN{1,i});
if OneNaN(1)==1 && range(rNaN{1,i})==0
    ConMito{1,i}(rNaN{i}( 1),:)=[]
else
for k=1:length(row1{1,i})
ConMito{1,i}((rNaNsort{1,i}(row1{1,i}(k), 1)-(k-1)), :)=[];
end
end



% % ConMito{1,i}(rNaN{i})=[];


% ConMito{1,i}(rN{i}( 1),:)=[];

 sContour_s1{1,i} = lssmooth([ConMito{1,i}(:,1),ConMito{1,i}(:,2)], 10 );
 sContour_s2{1,i} = lssmooth([ConMito{1,i}(:,3),ConMito{1,i}(:,4)], 10 );
 
 %calculate intensity density
 for r=1:(length(mesh{1,i})-1)
IntInt_av{1,i}(r,:)= (IntInt{1,i}(r,1)+IntInt{1,i}(r+1,1))/2;
IntDensity{1,i}(r,:)= IntInt_av{1,i}(r,:)/subsegmentArea{1,i}(r,:);
 end
 
MitoIntInt_av(1,i)=sum(IntInt_av{1,i}); %sum integrated intensity for full mito

 
%choose intensity density on mito slice which is most accurate
%measure integrated intensity over each mito to articially bleach slice
% MitoIntInt_av(i,1)=sum(IntInt_av{1,i});
% BleachFit




%find largest radii (top 5%) and their associated index
% [sort_radii{1,i}, indx{1,i}]=sort(radii{1,i}, 'descend');
% LenRad(1,i)=length(sort_radii{1,i});
% N_perc(1,i)=round(0.05*LenRad(1,i));
% maxRadii{1,i}=sort_radii{1,i}(1:N_perc(1,i));
% maxRadiiIndx{1,i}= indx{1,i}(1:N_perc(1,i));

%from these maxRadii, find the slice with the max integrated intensity
% [MaxIntInt(1,i),ChoiceS(1,i)] =max(IntInt_av{1,i}(maxRadiiIndx{1,i}));
% ChoiceSlice(1,i)=maxRadiiIndx{1,i}(ChoiceS(1,i));

%back-calculate area (thArea) and radius (thRad)for each section 
% thArea{1,i}=IntInt_av{1,i}(:,1)/IntDensity{1,i}(ChoiceSlice(1,i));

% % % % % % % % % % syms A B C x
% % % % % % % % % % A{1,i}=(1-RatioR{1,i}(r, :))*(1-RatioR{1,i}(r, :));
% % % % % % % % % % B{1,i}=meshSpacing*meshSpacing;
% % % % % % % % % % C{1,i}=thArea{1,i}(r, :)/(pi()*pi()*(RatioR{1,i}(r, :)+1)*(RatioR{1,i}(r, :)+1))
% % % % % % % % % % thRad{1,i}=solve(A*x^4+B*x^2+C);
% % % % % % % % % %  -(-(B - (B^2 - 4*A*C)^(1/2))/(2*A))^(1/2)
% % % % % % % % %   (-(B + (B^2 - 4*A*C)^(1/2))/(2*A))^(1/2)
% % % % % % % % %  -(-(B + (B^2 - 4*A*C)^(1/2))/(2*A))^(1/2)
% % % % % % % % %   (-(B - (B^2 - 4*A*C)^(1/2))/(2*A))^(1/2)
% % % % % % % % % OR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -(B - (B^2 - 4*A*C)^(1/2))/(2*A)
% -(B + (B^2 - 4*A*C)^(1/2))/(2*A)

% for rr=1:(length(mesh{1,i})-1)
% A{1,i}(rr,:)=(1-RatioR{1,i}(rr, :))*(1-RatioR{1,i}(rr, :));
% B{1,i}(rr,:)=meshSpacing*meshSpacing;
% C{1,i}(rr, :)=thArea{1,i}(rr, :)/(pi()*pi()*(RatioR{1,i}(rr, :)+1)*(RatioR{1,i}(rr, :)+1));
% 
% % Solution{1,i}(r,1)=-(-(B{1,i}(r,:) - (B{1,i}(r,:)^2 - 4*A{1,i}(r,:)*C{1,i}(r, :))^(1/2))/(2*A{1,i}(r,:)))^(1/2);
% % Solution{1,i}(r,2)=(-(B{1,i}(r,:) + (B{1,i}(r,:)^2 - 4*A{1,i}(r,:)*C{1,i}(r, :))^(1/2))/(2*A{1,i}(r,:)))^(1/2);
% % Solution{1,i}(r,3)=-(-(B{1,i}(r,:) + (B{1,i}(r,:)^2 - 4*A{1,i}(r,:)*C{1,i}(r, :))^(1/2))/(2*A{1,i}(r,:)))^(1/2);
% % Solution{1,i}(r,4)= (-(B{1,i}(r,:) - (B{1,i}(r,:)^2 - 4*A{1,i}(r,:)*C{1,i}(r, :))^(1/2))/(2*A{1,i}(r,:)))^(1/2);
% 
% Sol{1,i}(rr,1)=  (-B{1,i}(rr,:) - sqrt(B{1,i}(rr,:)*B{1,i}(rr,:) - 4*A{1,i}(rr,:)*C{1,i}(rr, :)))/(2*A{1,i}(rr,:));
% Sol{1,i}(rr,2)= (-B{1,i}(rr,:)+ sqrt(B{1,i}(rr,:)*B{1,i}(rr,:) - 4*A{1,i}(rr,:)*C{1,i}(rr, :)))/(2*A{1,i}(rr,:));
% 
% end

% for for r=1:(length(mesh{1,i})-1)
% syms A B C x
% A{1,i}(r,:)=(1-RatioR{1,i}(r,:))*(1-RatioR{1,i}(r,:));
% B{1,i}(r,:)=meshSpacing*meshSpacing;
% C{1,i}(r,:)=thArea{1,i}(r,:)/(pi()*pi()(RatioR{1,i}(r,:)+1)*(RatioR{1,i}(r,:)+1))
% thRad{1,i}(r,:)=solve(A*x^4+B*x^2+C);
%  end
 
% %choose intensity density on mito slice which is most accurate
% %find largest radii (top 5%) and their associated index
% [sort_radii{1,i}, indx{1,i}]=sort(radii{1,i}, 'descend');
% LenRad(1,i)=length(sort_radii{1,i});
% N_perc(1,i)=round(0.05*LenRad(1,i));
% maxRadii{1,i}=sort_radii{1,i}(1:N_perc(1,i));
% maxRadiiIndx{1,i}= indx{1,i}(1:N_perc(1,i));
% 
% %from these maxRadii, find the slice with the max integrated intensity
% [MaxIntInt(1,i),ChoiceS(1,i)] =max(IntInt_av{1,i}(maxRadiiIndx{1,i}));
% ChoiceSlice(1,i)=maxRadiiIndx{1,i}(ChoiceS(1,i));
% 
% %back-calculate area (thArea) for each section (assuming intensity density is
% %constant)
% thArea{1,i}(r,:)=




 
        end
 
    case 'No'
        exp=2;        
        %TO DO: call MitoMesh.m
end

%%intensity based contour (note that some necessary parameters are
%%calculated in the previous for loop)--- in progress

%remove infinite values from IntDensity cell, find index of min diameter
%for each frame
LocInf={};
ZeroLoc={};
MinLoc={};
MinVal={};

for ss=start_frame:nframe
         
LocInf{ss}=find(IntDensity{1,ss}==Inf);
IntDensity{1,ss}(LocInf{ss})=[];

ZeroLoc{ss}=find(radii{ss}==0);%avoid detection of diameter=0
radii{ss}(ZeroLoc{ss})=[];

LenRadii{ss}=length(radii{ss});
[MinVal{ss}, MinLoc{ss}]=min(radii{ss}(5:(LenRadii{ss}-5))); %to do: make user select a boundary over which to find the min

        
end

 
%manually select mito slice where intensity density is optimal
 imshow(cData(:,:,start_frame), []);
 hold on
 title('Choose (right click) 1 point on contour (uniform intensity, diffraction unlim, bright(ish))', 'fontsize', 12);
 plot(Contour_final{1,start_frame}(:,1),Contour_final{1,start_frame}(:,2), 'g-', 'linewidth', 1.2);
[SelectedPoint(:,1), SelectedPoint(:,2)]=getpts; %Col1 is X, Col2 Y, select first point only if user forgets to right click
 close;
 
%  imshow(cData(:,:,start_frame), []);
%  hold on
%  title('Choose (right click) 1 point on contour (at constriction site)', 'fontsize', 12);
%  plot(Contour_final{1,start_frame}(:,1),Contour_final{1,start_frame}(:,2), 'w-', 'linewidth', 1.2);
% [SelectedPointB(:,1), SelectedPointB(:,2)]=getpts; %Col1 is X, Col2 Y, select first point only if user forgets to right click
%  close;
 
rng=1; %range of mito over which bleaching is measured
for f=start_frame:nframe

  %find the slice (SelectSlice), closest to your first selection (this is
  %to measure a 'model' intensity density
 Xnorm=[Norm{1,start_frame}(:,1); Norm{1,start_frame}(:,3)];
 Ynorm=[Norm{1,start_frame}(:,2); Norm{1,start_frame}(:,4)];
 NormCat=[Xnorm Ynorm];
 
 IDX_norm=knnsearch(NormCat, SelectedPoint(1,:));
 
 del=IDX_norm-length(Norm{1,start_frame});
 
 %find the slice closest to your second choice (this is used to model
 %bleaching)
 
% IDX_normB=knnsearch(NormCat, SelectedPointB(1,:));
%  delB=IDX_normB-length(Norm{1,start_frame});
 
 
 
 if del<=0
   SelectSlice=Norm{1,start_frame}(IDX_norm,:); %coordinates of chosen slice (first frame)
    SelectIntDen=IntDensity{1,start_frame}(IDX_norm); %intensity density of chosen slice (first frame)
    SelectIntDenRg=sum(IntDensity{1,start_frame}(IDX_norm-rng:IDX_norm+rng)); %intensity density over range (first frame)
    for s=start_frame:nframe
        IntInt_Select{1,s}=IntInt_av{1,s}(IDX_norm-rng:IDX_norm+rng);
         
           IntDen_Select{1,s}=IntDensity{1,s}(IDX_norm-rng:IDX_norm+rng);
           IntDen_SelectCon{1,s}=IntDensity{1,s}(MinLoc{ss}-rng:MinLoc{ss}+rng);
    end
    
 MitoIntInt_Select=cellfun(@sum,IntInt_Select); %integrated intensity over selected region (the selected slice, including neighbours +/- rng)
%  MitoIntDen_Select=cellfun(@sum,IntDen_Select);%intensity density over selected region  (the selected slice, including neighbours +/- rng)
MitoIntDen_Select=cellfun(@mean,IntDen_Select);%mean of intensity density selected region of mito
%  TotMitoIntDen=cellfun(@sum, IntDensity); %sum of intensity density across full mito
TotMitoIntDen=cellfun(@mean, IntDensity); %mean of intensity density across full mito
 ConIntDen= cellfun(@sum, IntDen_SelectCon); %intensity density at constriction site
 id1=IDX_norm;
 
 else
      SelectSlice=Norm{1,start_frame}(del,:); %coordinates of slice
      SelectIntDen=IntDensity{1,start_frame}(del); %intensity density of chosen slice
      SelectIntDenRg=sum(IntDensity{1,start_frame}(del-rng:del+rng));
      for s=start_frame:nframe
        IntInt_Select{1,s}=IntInt_av{1,s}(del-rng:del+rng);

        IntDen_Select{1,s}=IntDensity{1,s}(del-rng:del+rng);
         IntDen_SelectCon{1,s}=IntDensity{1,s}(MinLoc{ss}-rng:MinLoc{ss}+rng);
    end
    MitoIntInt_Select=cellfun(@sum,IntInt_Select); %integrated intensity over selected region (the selected slice, including neighbours +/- rng)
%     MitoIntDen_Select=cellfun(@sum,IntDen_Select); %intensity density over selected region  (the selected slice, including neighbours +/- rng)
    
    MitoIntDen_Select=cellfun(@mean,IntDen_Select);
    TotMitoIntDen=cellfun(@mean, IntDensity);
    
%     TotMitoIntDen=cellfun(@sum, IntDensity);%intensity density across full mito
    ConIntDen= cellfun(@sum, IntDen_SelectCon); %intensity density at constriction site
    id1=del
 end
 
% % %   if delB<=0
% % %     SelectSliceB=Norm{1,start_frame}(IDX_normB,:); %coordinates of slice (first frame)
% % %     SelectIntDenB=IntDensity{1,start_frame}(IDX_normB); %intensity density of chosen slice (first frame)
% % %   for s=start_frame:nframe
% % %         IntInt_SelectB{1,s}=IntInt_av{1,s}(IDX_normB-rng:IDX_normB+rng);
% % %          
% % %         IntDen_SelectB{1,s}=IntDensity{1,s}(IDX_normB-rng:IDX_normB+rng);
% % %         
% % %   end
% % %      MitoIntDen_SelectB=cellfun(@sum,IntDen_SelectB);%intensity density over selected region
% % %      id2=IDX_normB;
% % %   else
% % %       SelectSliceB=Norm{1,start_frame}(delB,:); %coordinates of slice (first frame)
% % %     SelectIntDenB=IntDensity{1,start_frame}(delB); %intensity density of chosen slice (first frame)
% % %   for s=start_frame:nframe
% % %         IntInt_SelectB{1,s}=IntInt_av{1,s}(delB-rng:delB+rng);
% % %          
% % %         IntDen_SelectB{1,s}=IntDensity{1,s}(delB-rng:delB+rng);
% % %         
% % %   end
% % %      MitoIntDen_SelectB=cellfun(@sum,IntDen_SelectB);%intensity density over selected region
% % %      id2=delB
% % %   end
% % %   
% % %  


 %measure/apply bleaching 
 
% MitoIntInt_Select(1,f)=sum(IntInt_av{1,f}(del-rng:del+rng, :)); %sum integrated intensity for part of mito
BleachFit=polyfit((start_frame:1:nframe)', MitoIntDen_Select(start_frame:end)', 1); % bleaching of selected int density region
% BleachFit=polyfit((start_frame:1:nframe)',  TotMitoIntDen(start_frame:end)', 1); %bleaching of total mito intensity density
% BleachFit=polyfit((start_frame:1:nframe)',  ConIntDen(start_frame:end)', 1); %bleaching around constriction site
RateOfChange=BleachFit(1,1); %slope of bleach curve for full movie (IntDen/frame)
intercept=BleachFit(1,2);
 
%   BleachedIntDen(1, start_frame)=SelectIntDen;
%   BleachedIntDen(1,f+1)= BleachedIntDen(1,f)+RateOfChange; %bleached slices, not that change in frames is always 1
%  
 % BleachedIntDen(1, start_frame)=MitoIntDen_Select(start_frame);
%   Inter=MitoIntDen_Select(start_frame)-(start_frame*RateOfChange);
% Inter=SelectIntDen;
Inter=SelectIntDen-RateOfChange*(start_frame);
BleachedIntDen(1,f)=(RateOfChange*f)+Inter;
%  BleachedIntDen(1,f)=(RateOfChange*f)+intercept;
 
 %backcalculate the area (thArea) of each slice 
%  thArea{1,f}=IntInt_av{1,f}(:,1)/BleachedIntDen(1,start_frame); %this does not take bleaching into account
 thArea{1,f}=IntInt_av{1,f}(:,1)/BleachedIntDen(1,f); %take into account bleaching
 %from the backcalculated areas, find corresponding radii
 for rr=1:(length(mesh{1,f})-1)
A{1,f}(rr,:)=(1-RatioR{1,f}(rr, :))*(1-RatioR{1,f}(rr, :));
B{1,f}(rr,:)=meshSpacing*meshSpacing;
C{1,f}(rr, :)=(thArea{1,f}(rr, :))^2/(pi()*pi()*(RatioR{1,f}(rr, :)+1)*(RatioR{1,f}(rr, :)+1));

Sol{1,f}(rr,1)=  (-B{1,f}(rr,:) - sqrt(B{1,f}(rr,:)*B{1,f}(rr,:) - 4*A{1,f}(rr,:)*C{1,f}(rr, :)))/(2*A{1,f}(rr,:));
Sol{1,f}(rr,2)= (-B{1,f}(rr,:) + sqrt(B{1,f}(rr,:)*B{1,f}(rr,:) - 4*A{1,f}(rr,:)*C{1,f}(rr, :)))/(2*A{1,f}(rr,:)); %take magnitude of real part of this solution

end
 
end
   
%calculate the curvature of entire contour for every frame (performed on
%FWHM contour below)

for jj=start_frame:nframe
    %side1
   ContnextX1{:,jj}= sContour_s1{:,jj}(:,1);
   ContnextY1{:,jj}= sContour_s1{:,jj}(:,2);
   v1{jj}=[ContnextX1{:,jj}, ContnextY1{:,jj}];
   Lines1{jj}=[(1:(size(v1{1,jj},1)-2))' (3:size(v1{1,jj},1))'];% the curvature at a single point uses 3 neighbouring points (7 total: itself and 3 on each side)
    
  [k_01{jj}]=LineCurvature2D(v1{jj}, Lines1{jj});
  [k1{jj}]=[k_01{jj}./px];
  
   %side2
   ContnextX2{:,jj}= sContour_s2{:,jj}(:,1);
   ContnextY2{:,jj}= sContour_s2{:,jj}(:,2);
   v2{jj}=[ContnextX2{:,jj}, ContnextY2{:,jj}];
   Lines2{jj}=[(1:(size(v2{1,jj},1)-2))' (3:size(v2{1,jj},1))'];% the curvature at a single point uses 3 neighbouring points (7 total: itself and 3 on each side)
    
  [k_02{jj}]=LineCurvature2D(v2{jj}, Lines2{jj});
  [k2{jj}]=[k_02{jj}./px];
  
end

%visually inspect contour over SIM data
for ii=start_frame:nVisIn 
   h=figure;
   axis equal
   imshow(cData(:,:,ii), []);
   hold on;
   plot(sContour_s1{1,ii}(:,1), sContour_s1{1,ii}(:,2), 'y-', 'linewidth', 1.8);
   plot(sContour_s2{1,ii}(:,1), sContour_s2{1,ii}(:,2), 'y-', 'linewidth', 1.8);
   plot(Contour_final{1,ii}(:,1),Contour_final{1,ii}(:,2), 'm-', 'linewidth', 1.2);
   plot(mesh{1,ii}(id1, 1:2:3), mesh{1,ii}(id1,2:2:4), 'g', 'linewidth', 1.5);
%     plot(mesh{1,ii}(id2, 1:2:3), mesh{1,ii}(id2,2:2:4), 'w', 'linewidth', 1.5);
%    plot(ConMito{1,ii}(:,3),ConMito{1,ii}(:,4), 'g-', 'linewidth', 1.2);
%      plot(ConMito{1,ii}(:,1),ConMito{1,ii}(:,2), 'g-', 'linewidth', 1.2);
end


        



%%%TO DO: 
%1-redefine ChoiceSlice
%2-choose 'best' solution based on back-calculated radius
%3-shift contour based on 2





end

