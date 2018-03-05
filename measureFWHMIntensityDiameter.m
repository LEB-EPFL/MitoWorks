function [ContourFWHMside1 ContourFWHMside2 DiametersFWHM DiametersIntensity, ConMito, dist, Int, fitG] = measurefitFWHMIntensityDiameter(OriginalImage,pos,subNormals,mesh,meshSpacing, ContnextXsnake, ContnextYsnake,StartFrame, EndFrame)
%generateFinalContours produces a final mitochondria contour and 2 contours for
%comparison based on the SIM image

%   Detailed explanation goes here
%script that generates a mitochondrial contour, calculates
%envelope curvature along contour and measures distances across mito
%according to the intensity histograms




for f=StartFrame:EndFrame
    Image{f} = OriginalImage{f};
  % Image{f}= im2uint16(Image{f});
     
    I(:,:,f)=Image{f};
    %Icrop(:,:,f)=imcrop(I(:,:,f),(pos./3)); 
    
end

   

    for i=StartFrame:EndFrame
    for j=1:length(subNormals{i})
        
   subNormals{:,i}(all(subNormals{:,i}==0,2), :)=[];  
   
   TotDist{:,i}=sqrt((subNormals{:,i}(:,1)-subNormals{:,i}(:,3)).*(subNormals{:,i}(:,1)-subNormals{:,i}(:,3))+(subNormals{:,i}(:,2)-subNormals{:,i}(:,4)).*(subNormals{:,i}(:,2)-subNormals{:,i}(:,4)));
    end
    end
    
%     subNormals=cellfun(@(x) x*px,subNormals1,'un',0);
% h = waitbar(0,'Measuring FWHM diameter values...');

 for i=StartFrame:EndFrame
%              waitbar((i-StartFrame)/(EndFrame-StartFrame));

    for j=1:length(subNormals{i})
        
[cx{j,i}, cy{j,i}, c{j,i}, xi{j,i}, yi{j,i}]=improfile(I(:,:,i), [subNormals{1,i}(j,1), subNormals{1,i}(j,3)], [subNormals{1,i}(j,2), subNormals{1,i}(j,4)], 'bilinear' );
    end
 end
 ValleyL2=cell(size(c,1),EndFrame);
IDleft=cell(size(c,1),EndFrame);
idLeftBound=cell(size(c,1),EndFrame);
ValleyR2=cell(size(c,1),EndFrame);
IDright=cell(size(c,1),EndFrame);
idRightBound=cell(size(c,1),EndFrame);
%cx, cy are the spatial coordinates in pixels across the profile, c is the
%intensity of the profile, xi and yi are the inital points
 for i=StartFrame:EndFrame
         for j=1:length(subNormals{i})
%subtract background from intensity c{j,i} 
MinI{j,i}=min(c{j,i});
Int{j,i}=c{j,i}-MinI{j,i}; %background subtracted profile
Int{j,i}(isnan(Int{j,i}(:,1)),:)=[];
% Int{j,i}=lssmooth(Int{j,i},4);

%transform real space coordinates to distance along line and compute
%integrated intensity for each profile
Dx{j,i}=cx{j,i}-xi{j,i}(1,1);
Dy{j,i}=cy{j,i}-yi{j,i}(1,1);
dist{j,i}= sqrt((Dx{j,i}.*Dx{j,i})+(Dy{j,i}.*Dy{j,i})); %this is the x distance along x-axis
dist{j,i}(isnan(dist{j,i}(:,1)),:)=[];

if size(Int{j,i},1)==size(dist{j,i},1);
IntInt{1,i}(j,:)=trapz(dist{j,i}, Int{j,i}); %integrated intensity for each profile
end

%calculate truncated cone area
radii{1,i}=sqrt((mesh{1,i}(:,1)-mesh{1,i}(:,3)).^2+(mesh{1,i}(:,2)-mesh{1,i}(:,4)).^2)/2;
for r=1:(length(mesh{1,i})-1)
subsegmentArea{1,i}(r,:)=(pi()*(radii{1,i}(r,1)+radii{1,i}(r+1,1))*sqrt((radii{1,i}(r,1)-radii{1,i}(r+1,1))^2+(meshSpacing)^2));
subsegmentVol{1,i}(r, :)=(pi()*(1/3)*(meshSpacing)*((radii{1,i}(r,1))^2+(radii{1,i}(r+1,1))^2)+radii{1,i}(r+1,1)*radii{1,i}(r,1));
RatioR{1,i}(r, :)=radii{1,i}(r+1,1)/radii{1,i}(r,1);

% IntInt_av{1,i}(r,:)= (IntInt{1,i}(r,1)+IntInt{1,i}(r+1,1))/2;
% Int_Dis{1,i}(r,:)= IntInt_av{1,i}(r,:)./IntInt_av{1,i}(r,:);
end


% MitoIntInt_av(1,i)=sum(IntInt_av{1,i}); %sum integrated intensity for full mito
% %perform gaussian fit
% t=~isnan(Int{j,i})& ~isnan(dist{j,i});
if size(Int{j,i},1)==size(dist{j,i},1);
    
[pks{j,i},locs{j,i},width{j,i}, prominence{j,i}] = findpeaks(Int{j,i},dist{j,i}); %pks and locs are the Y and X values of the peaks
[MaxPks{j,i}, Maxrow{1,i}(:,j)]=max(pks{j,i}); %MaxPks is the max peak and Maxrow is the row at which its found in the pks matrix
IDpk{j,i}=find(Int{j,i}==MaxPks{j,i}); %IDpk is the index of the peak in the dist and intensity matrices
[valley{j,i},locsVal{j,i}]=findpeaks((-Int{j,i}), dist{j,i}); %valley and locsVal are the -Y and X values of the valleys
if isempty(locsVal{j,i})==0
% for p=1:size(locsVal{j,i},1)
% locsVal{j,i}<dist{j,i}(IDpk{j,i})



ValleyL1{j,i}=[locsVal{j,i}(locsVal{j,i}<dist{j,i}(IDpk{j,i})),-valley{j,i}(locsVal{j,i}<dist{j,i}(IDpk{j,i}))]; %define right valley
if isempty(ValleyL1{j,i})==0
ValleyL2{j,i}=[ValleyL1{j,i}(ValleyL1{j,i}(:,2)<=0.7*abs(MaxPks{j,i}),1),ValleyL1{j,i}(ValleyL1{j,i}(:,2)<=0.7*abs(MaxPks{j,i}),2)]; %filter left valley values that are within 30% of max peak
end

% ValleyL{j,i}=[locsVal{j,i}(abs(valley{j,i})<=0.7*abs(MaxPks{j,i})),-valley{j,i}(abs(valley{j,i})<=0.7*abs(MaxPks{j,i}))]; %filter left valley values that are within 30% of max peak
% ValleyL{j,i}=[locsVal{j,i}(locsVal{j,i}<dist{j,i}(IDpk{j,i})),-valley{j,i}(locsVal{j,i}<dist{j,i}(IDpk{j,i}))]; %define left valley as all valleys to left of peak
% ValleyL{j,i}=[locsVal{j,i}(locsVal{j,i}<(dist{j,i}(IDpk{j,i})-0.5*(width{j,i}(Maxrow{1,i}(:,j))))),-valley{j,i}(locsVal{j,i}<(dist{j,i}(IDpk{j,i})-0.5*(width{j,i}(Maxrow{1,i}(:,j)))))];
%find closest left valley point to maximum peak

if isempty(ValleyL2{j,i})==0
IDleft{j,i}=knnsearch(ValleyL2{j,i},[dist{j,i}(IDpk{j,i}),Int{j,i}(IDpk{j,i})],'K',1 ); %this is the row index of the valley in the ValleyL matrix
if isempty(IDleft{j,i})==0
idLeftBound{j,i}=find(dist{j,i}==ValleyL2{j,i}(IDleft{j,i},1)); %this is the ID of the left valley in the Int, dist matrices
end
%set everything to the left of the left valley to the Intensity value of the left valley
% Int{j,i}(dist{j,i}<ValleyL{j,i}(IDleft{j,i}),1)=ValleyL{j,i}(IDleft{j,i},2);
% MinI2{j,i}=min(Int{j,i});
% Int{j,i}=Int{j,i}-MinI2{j,i};

if isempty(idLeftBound{j,i})==0
Int{j,i}(dist{j,i}<dist{j,i}(idLeftBound{j,i}))=[];
dist{j,i}(dist{j,i}<dist{j,i}(idLeftBound{j,i}))=[];

end
end
%redefine the location of max peak in dist and Int matrices since these
%changed once Int and dist changed
[pksR{j,i},locsR{j,i},widthR{j,i}, prominenceR{j,i}] = findpeaks(Int{j,i},dist{j,i}); %pks and locs are the Y and X values of the peaks
[MaxPksR{j,i}, MaxrowR{1,i}(:,j)]=max(pksR{j,i}); %MaxPks is the max peak and Maxrow is the row at which its found in the pks matrix
IDpkR{j,i}=find(Int{j,i}==MaxPksR{j,i}); %IDpk is the index of the peak in the dist and intensity matrices



% if locsVal{j,i}>dist{j,i}(IDpk{j,i})
ValleyR1{j,i}=[locsVal{j,i}(locsVal{j,i}>dist{j,i}(IDpkR{j,i})),-valley{j,i}(locsVal{j,i}>dist{j,i}(IDpkR{j,i}))]; %define right valley


if isempty(ValleyR1{j,i})==0
ValleyR2{j,i}=[ValleyR1{j,i}(ValleyR1{j,i}(:,2)<=0.7*abs(MaxPksR{j,i}),1),ValleyR1{j,i}(ValleyR1{j,i}(:,2)<=0.7*abs(MaxPksR{j,i}),2)]; %filter left valley values that are within 30% of max peak
end
% ValleyR{j,i}=[locsVal{j,i}(locsVal{j,i}>(dist{j,i}(IDpk{j,i})+0.5*(width{j,i}(Maxrow{1,i}(:,j))))),-valley{j,i}(locsVal{j,i}>(dist{j,i}(IDpk{j,i})+0.5*(width{j,i}(Maxrow{1,i}(:,j)))))];
if isempty(ValleyR2{j,i})==0
IDright{j,i}=knnsearch(ValleyR2{j,i},[dist{j,i}(IDpkR{j,i}),Int{j,i}(IDpkR{j,i})], 'K',1);
if isempty(IDright{j,i})==0
idRightBound{j,i}=find(dist{j,i}==ValleyR2{j,i}(IDright{j,i},1));
end
% Int{j,i}(dist{j,i}>ValleyR{j,i}(IDright{j,i}),1)=ValleyR{j,i}(IDright{j,i},2);
if isempty(idRightBound{j,i})==0
 Int{j,i}(dist{j,i}>dist{j,i}(idRightBound{j,i}))=[];
dist{j,i}(dist{j,i}>dist{j,i}(idRightBound{j,i}))=[];
end
end
% MinI2{j,i}=min(Int{j,i});
% Int{j,i}=Int{j,i}-MinI2{j,i};
end

%find higher min and subtract from intensity profile
% MinL{j,i}=min(Int{j,i}(dist{j,i}<dist{j,i}(IDpk{j,i})));

% MinR{j,i}=min(Int{j,i}(dist{j,i}>dist{j,i}(IDpk{j,i})));
% AbsMin{j,i}=max([MinL{j,i};MinR{j,i}]);




% options{j,i}=fitoptions('gauss1');
%    options{j,i}.Lower=[0, dist{j,i}(IDpk{j,i})-0.1, 0];
%    options{j,i}.Upper=[Inf,dist{j,i}(IDpk{j,i})+0.1, Inf ];
%     fitG{j,i}=fit(dist{j,i}, (Int{j,i}-AbsMin{j,i}), 'gauss1');


    fitG{j,i}=fit(dist{j,i}, (Int{j,i}-min(Int{j,i})), 'gauss1');
    Param{j,i}=coeffvalues(fitG{j,i});
    FWHM{j,i}=2.354*(1/sqrt(2))*Param{j,i}(3);
    CentreFit{j,i}=Param{j,i}(2);

end
 
    
    end
 end
 

% %find the FWHM points of intersection with x-axis
% % % [width{j,i}, xLeft{j,i}, xRight{j,i},yLeft{j,i},yRight{j,i}]= getFWHM(dist{j,i},Int{j,i});
% [width{1,i}(j, :), xLeft{1,i}(j, :), xRight{1,i}(j,:),yLeft{1,i}(j, :),yRight{1,i}(j,:)]= getFWHM(dist{j,i},c{j,i});
for i=StartFrame:EndFrame
for j=1:length(CentreFit(:,i))
    
if isempty(CentreFit{j,i})==0
xLeft{j,i}=CentreFit{j,i}-0.5*FWHM{j,i};
xRight{j,i}=CentreFit{j,i}+0.5*FWHM{j,i};



% % % DiametersFWHM{1,i}(j,:)=width{j,i};
DiametersFWHM{1,i}(j,:)=FWHM{j,i};

% %transform back into real space
% IDalign{j,i}=knnsearch([xi{j,i},yi{j,i}],[subNormals{1,i}(j,1),subNormals{1,i}(j,2)],'K',1);
% if IDalign{j,i}==1
% IDother{j,i}=2;
% else
%     IDother{j,i}=1;   
% end
% delta_x{j,i}=xi{j,i}(IDother{j,i})-xi{j,i}(IDalign{j,i});
delta_x{j,i}=xi{j,i}(2)-xi{j,i}(1);
delta_y{j,i}=yi{j,i}(2)-yi{j,i}(1);
% delta_y{j,i}=yi{j,i}(IDother{j,i})-yi{j,i}(IDalign{j,i});
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
 end
    end
end
end

 
%  close(h);
%  DiametersFWHM=width{j,i};
    


 %remove NaN values 
%  for i=StartFrame: EndFrame
% FindNaN= cellfun(@(x) (isnan(x)),ConMito, 'UniformOutput', false);
% [rNaN{i}, cNaN{i}]=find(FindNaN{1,i}==1);
% rNaNsort{1,i}=sort(rNaN{1,i});
% rNshift{1,i}=circshift(rNaNsort{1,i},1);
% [row1{1,i}, col1{1,i}]=find(abs(rNshift{1,i}(:,1)- rNaNsort{1,i}(:,1))==1);
% 
% OneNaN=any(FindNaN{1,i});
% if OneNaN(1)==1 && range(rNaN{1,i})==0
%     ConMito{1,i}(rNaN{i}( 1),:)=[];
% else
% for k=1:length(row1{1,i})
% ConMito{1,i}((rNaNsort{1,i}(row1{1,i}(k), 1)-(k-1)), :)=[];
% end
% end
% 
% OneNaN=any(FindNaN{1,i});
% if OneNaN(1)==1 && range(rNaN{1,i})==0
%     ConMito{1,i}(rNaN{i}( 1),:)=[];
% end








% % ConMito{1,i}(rNaN{i})=[];


% ConMito{1,i}(rN{i}( 1),:)=[];
for i=StartFrame:EndFrame
 ContourFWHMside1{1,i} = lssmooth([ConMito{1,i}(:,1),ConMito{1,i}(:,2)], 6);
 ContourFWHMside2{1,i} = lssmooth([ConMito{1,i}(:,3),ConMito{1,i}(:,4)], 6);
end
 
 %calculate intensity density
 for i=StartFrame:EndFrame
 for r=1:(length(ContourFWHMside1{1,i})-1) %%warning: matrices from these loops, might be shifted from the mesh indices
%      if isempty(nonzeros(mesh{1,i}(1)))==1
IntInt_av{1,i}(r,:)= (IntInt{1,i}(r,1)+IntInt{1,i}(r+1,1))/2;
IntDensity{1,i}(r,:)= IntInt_av{1,i}(r,:)/subsegmentArea{1,i}(r,:);
IntDensityVol{1,i}(r,:)= IntInt_av{1,i}(r,:)/subsegmentVol{1,i}(r,:);
     
 end
 
MitoIntInt_av(1,i)=sum(IntInt_av{1,i}); %sum integrated intensity for full mito
% end
 end


 
        
 
%     case 'No'
%         exp=2;        
        %TO DO: call MitoMesh.m
%end

%%intensity based measurement (note that some necessary parameters are
%%calculated in the previous for loop)--- in progress

%remove infinite values from IntDensity cell, find index of min diameter
%for each frame
LocInf={};
ZeroLoc={};
MinLoc={};
MinVal={};

for ss=StartFrame:EndFrame
         
LocInf{ss}=find(IntDensity{1,ss}==Inf);
IntDensity{1,ss}(LocInf{ss})=[];
IntDensityVol{1,ss}(LocInf{ss})=[];

ZeroLoc{ss}=find(radii{ss}==0);%avoid detection of diameter=0
radii{ss}(ZeroLoc{ss})=[];

LenRadii{ss}=length(radii{ss});
[MinVal{ss}, MinLoc{ss}]=min(radii{ss}(5:(LenRadii{ss}-5))); %to do: make user select a boundary over which to find the min

        
end

 
%manually select mito slice where intensity density is optimal
%  imshow(I(:,:,StartFrame), []);
%  hold on
%  title('Choose (right click) 1 point on contour (uniform intensity, diffraction unlim, bright(ish))', 'fontsize', 12);
%  plot(ContnextXsnake{1,StartFrame}(:,1),ContnextYsnake{1,StartFrame}(:,1), 'g-', 'linewidth', 1.2);
% [SelectedPoint(:,1), SelectedPoint(:,2)]=getpts; %Col1 is X, Col2 Y, select first point only if user forgets to right click
%  close;
 

 
% rng=1; %range of mito over which model Intensity density is defined
for f=StartFrame:EndFrame
% 
%   %find the slice (SelectSlice), closest to your first selection (this is
%   %to measure a 'model' intensity density
%  XsubNormals=[subNormals{1,StartFrame}(:,1); subNormals{1,StartFrame}(:,3)];
%  YsubNormals=[subNormals{1,StartFrame}(:,2); subNormals{1,StartFrame}(:,4)];
%  subNormalsCat=[XsubNormals YsubNormals];
%  
%  IDX_subNormals=knnsearch(subNormalsCat, SelectedPoint(1,:));
%  
%  del=IDX_subNormals-length(subNormals{1,StartFrame});
%  
 %find the slice closest to your second choice (this is used to model
 %bleaching)
 
% IDX_subNormalsB=knnsearch(subNormalsCat, SelectedPointB(1,:));
%  delB=IDX_subNormalsB-length(subNormals{1,StartFrame});
 
 
 
%  if del<=0
%    SelectSlice=subNormals{1,StartFrame}(IDX_subNormals,:); %coordinates of chosen slice (first frame)
%     SelectIntDen=IntDensity{1,StartFrame}(IDX_subNormals); %intensity density of chosen slice (first frame)
%     SelectIntDenVol=IntDensityVol{1,StartFrame}(IDX_subNormals); %intensity density of chosen slice (first frame)
% 
%     SelectIntDenRg=sum(IntDensity{1,StartFrame}(IDX_subNormals-rng:IDX_subNormals+rng)); %intensity density over range (first frame)
%     SelectIntDenVolRg=sum(IntDensityVol{1,StartFrame}(IDX_subNormals-rng:IDX_subNormals+rng)); %intensity density over range (first frame)
% 
%     for s=StartFrame:EndFrame
%         IntInt_Select{1,s}=IntInt_av{1,s}(IDX_subNormals-rng:IDX_subNormals+rng);
%          
%            IntDen_Select{1,s}=IntDensity{1,s}(IDX_subNormals-rng:IDX_subNormals+rng);
%            IntDenVol_Select{1,s}=IntDensityVol{1,s}(IDX_subNormals-rng:IDX_subNormals+rng);
%             
%            %selected intensity density around constriction  site
%            IntDen_SelectCon{1,s}=IntDensity{1,s}(MinLoc{ss}-rng:MinLoc{ss}+rng); %int density based on area
%            IntDenVol_SelectCon{1,s}=IntDensityVol{1,s}(MinLoc{ss}-rng:MinLoc{ss}+rng); %int density based on vol
% 
%     end
    
%  MitoIntInt_Select=cellfun(@sum,IntInt_Select); %integrated intensity over selected region (the selected slice, including neighbours +/- rng)
% %  MitoIntDen_Select=cellfun(@sum,IntDen_Select);%intensity density over selected region  (the selected slice, including neighbours +/- rng)
% MitoIntDen_Select=cellfun(@mean,IntDen_Select);%mean of intensity density selected region of mito
% MitoIntDenVol_Select=cellfun(@mean,IntDenVol_Select);%mean of intensity density (based on volume) at selected region of mito

%  TotMitoIntDen=cellfun(@sum, IntDensity); %sum of intensity density across full mito
TotMitoIntDen=cellfun(@mean, IntDensity); %mean of intensity density across full mito
TotMitoIntDenVol=cellfun(@mean, IntDensityVol); %mean of intensity density across full mito

%  ConIntDen= cellfun(@sum, IntDen_SelectCon); %intensity density based on Area at constriction site
%   ConIntDenVol= cellfun(@sum, IntDenVol_SelectCon); %intensity density based on volume at constriction site

%  id1=IDX_subNormals;
%  
%  else
%       SelectSlice=subNormals{1,StartFrame}(del,:); %coordinates of slice
%       SelectIntDen=IntDensity{1,StartFrame}(del); %intensity density of chosen slice
%       SelectIntDenVol=IntDensityVol{1,StartFrame}(del); %intensity density based on volume of chosen slice
% 
%       SelectIntDenRg=sum(IntDensity{1,StartFrame}(del-rng:del+rng));
%       SelectIntDenVolRg=sum(IntDensityVol{1,StartFrame}(del-rng:del+rng));
%       for s=StartFrame:EndFrame
%         IntInt_Select{1,s}=IntInt_av{1,s}(del-rng:del+rng);
% 
%         IntDen_Select{1,s}=IntDensity{1,s}(del-rng:del+rng);
%         IntDenVol_Select{1,s}=IntDensityVol{1,s}(del-rng:del+rng);
%          IntDen_SelectCon{1,s}=IntDensity{1,s}(MinLoc{ss}-rng:MinLoc{ss}+rng);
%          IntDenVol_SelectCon{1,s}=IntDensityVol{1,s}(MinLoc{ss}-rng:MinLoc{ss}+rng);
%     end
%     MitoIntInt_Select=cellfun(@sum,IntInt_Select); %integrated intensity over selected region (the selected slice, including neighbours +/- rng)
%     MitoIntDen_Select=cellfun(@sum,IntDen_Select); %intensity density over selected region  (the selected slice, including neighbours +/- rng)
    
%     MitoIntDen_Select=cellfun(@mean,IntDen_Select);
%      MitoIntDenVol_Select=cellfun(@mean,IntDenVol_Select);
    TotMitoIntDen=cellfun(@mean, IntDensity);
    TotMitoIntDenVol=cellfun(@mean, IntDensityVol);

    
%     TotMitoIntDen=cellfun(@sum, IntDensity);%intensity density across full mito
%     ConIntDen= cellfun(@sum, IntDen_SelectCon); %intensity density at constriction site
%        ConIntDenVol= cellfun(@sum, IntDenVol_SelectCon); %intensity density at constriction site
%     id1=del;
%  end
 %%%%%you are here with edits to vol
 %measure/apply bleaching 
 
% MitoIntInt_Select(1,f)=sum(IntInt_av{1,f}(del-rng:del+rng, :)); %sum integrated intensity for part of mito
% BleachFit=polyfit((StartFrame:1:EndFrame)', MitoIntDen_Select(StartFrame:end)', 1); % bleaching of selected int density region
% BleachFit=polyfit((StartFrame:1:EndFrame)',  TotMitoIntDen(StartFrame:end)', 1); %bleaching of total mito intensity density based on area
BleachFit=polyfit((StartFrame:1:EndFrame)',  TotMitoIntDenVol(StartFrame:end)', 1); %bleaching of total mito intensity density based on volume

% BleachFit=polyfit((StartFrame:1:EndFrame)',  ConIntDen(StartFrame:end)', 1); %bleaching around constriction site
RateOfChange=BleachFit(1,1); %slope of bleach curve for full movie (IntDen/frame)
intercept=BleachFit(1,2);
 
%   BleachedIntDen(1, StartFrame)=SelectIntDen;
%   BleachedIntDen(1,f+1)= BleachedIntDen(1,f)+RateOfChange; %bleached slices, not that change in frames is always 1
%  
 % BleachedIntDen(1, StartFrame)=MitoIntDen_Select(StartFrame);
%   Inter=MitoIntDen_Select(StartFrame)-(StartFrame*RateOfChange);
% Inter=SelectIntDen;
% Inter=SelectIntDen-RateOfChange*(StartFrame);
% BleachedIntDen(1,f)=(RateOfChange*f)+Inter;
 BleachedIntDen(1,f)=(RateOfChange*f)+intercept;
 
 %backcalculate the area (thArea) of each slice 
%  thArea{1,f}=IntInt_av{1,f}(:,1)/BleachedIntDen(1,StartFrame); %this does not take bleaching into account
%  thArea{1,f}=IntInt_av{1,f}(:,1)/BleachedIntDen(1,f); %take into account
%  bleaching to calculate the area based on intensity density
 thVol{1,f}=IntInt_av{1,f}(:,1)/BleachedIntDen(1,f);%take into account bleaching
 %from the backcalculated areas, find corresponding radii
% % %  for rr=1:(length(mesh{1,f})-1)
% % % A{1,f}(rr,:)=(1-RatioR{1,f}(rr, :))*(1-RatioR{1,f}(rr, :));
% % % B{1,f}(rr,:)=meshSpacing*meshSpacing;
% % % C{1,f}(rr, :)=(thArea{1,f}(rr, :))^2/(pi()*pi()*(RatioR{1,f}(rr, :)+1)*(RatioR{1,f}(rr, :)+1));
% % % 
% % % Sol{1,f}(rr,1)=  (-B{1,f}(rr,:) - sqrt(B{1,f}(rr,:)*B{1,f}(rr,:) - 4*A{1,f}(rr,:)*C{1,f}(rr, :)))/(2*A{1,f}(rr,:));
% % % Sol{1,f}(rr,2)= (-B{1,f}(rr,:) + sqrt(B{1,f}(rr,:)*B{1,f}(rr,:) - 4*A{1,f}(rr,:)*C{1,f}(rr, :)))/(2*A{1,f}(rr,:)); %take square root of the magnitude of real part this solution
% % % DiametersIntensity{1,f}=sqrt(abs(Sol{1,f}(:,2)))
% % % end

%from the backcalculated volumes, find corresponding radii
for rr=1:(length(ContourFWHMside1{1,f})-1)
     Radius_Sol{1,f}(rr, :)=sqrt(thVol{1,f}(rr, :)*3*(1/meshSpacing)*(1/pi())*(1+(RatioR{1,f}(rr,:))^2+RatioR{1,f}(rr,:))^(-1));
    DiametersIntensity{1,f}(rr,:)=2* Radius_Sol{1,f}(rr, :);
end
end
   
%calculate the curvature of entire contour for every frame (performed on
%FWHM contour below)

% for jj=StartFrame:EndFrame
%     %side1
%    ContnextX1{:,jj}= sContour_s1{:,jj}(:,1);
%    ContnextY1{:,jj}= sContour_s1{:,jj}(:,2);
%    v1{jj}=[ContnextX1{:,jj}, ContnextY1{:,jj}];
%    Lines1{jj}=[(1:(size(v1{1,jj},1)-2))' (3:size(v1{1,jj},1))'];% the curvature at a single point uses 3 neighbouring points (7 total: itself and 3 on each side)
%     
%   [k_01{jj}]=LineCurvature2D(v1{jj}, Lines1{jj});
%   [k1{jj}]=[k_01{jj}./px];
%   
%    %side2
%    ContnextX2{:,jj}= sContour_s2{:,jj}(:,1);
%    ContnextY2{:,jj}= sContour_s2{:,jj}(:,2);
%    v2{jj}=[ContnextX2{:,jj}, ContnextY2{:,jj}];
%    Lines2{jj}=[(1:(size(v2{1,jj},1)-2))' (3:size(v2{1,jj},1))'];% the curvature at a single point uses 3 neighbouring points (7 total: itself and 3 on each side)
%     
%   [k_02{jj}]=LineCurvature2D(v2{jj}, Lines2{jj});
%   [k2{jj}]=[k_02{jj}./px];
%   
% end

%visually inspect contour over SIM data
for ii=StartFrame:EndFrame
 h=figure;
   axis equal
   imshow(I(:,:,ii), []);
   hold on;
   plot(ContourFWHMside1{1,ii}(:,1), ContourFWHMside1{1,ii}(:,2), 'g-', 'linewidth', 1.4);
   plot(ContourFWHMside2{1,ii}(:,1), ContourFWHMside2{1,ii}(:,2), 'g-', 'linewidth', 1.4);
   plot(ContnextXsnake{1,ii}(:,1),ContnextYsnake{1,ii}(:,1), 'm-', 'linewidth', 1.2);
%    plot(mesh{1,ii}(id1, 1:2:3), mesh{1,ii}(id1,2:2:4), 'g', 'linewidth', 1.5);
   drawnow;
%     plot(mesh{1,ii}(id2, 1:2:3), mesh{1,ii}(id2,2:2:4), 'w', 'linewidth', 1.5);
%    plot(ConMito{1,ii}(:,3),ConMito{1,ii}(:,4), 'g-', 'linewidth', 1.2);
%      plot(ConMito{1,ii}(:,1),ConMito{1,ii}(:,2), 'g-', 'linewidth', 1.2);
end


        



%%%TO DO: 
%1-redefine ChoiceSlice
%2-choose 'best' solution based on back-calculated radius
%3-shift contour based on 2





% end

