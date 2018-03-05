function [ContourFWHMside1 ContourFWHMside2 DiametersFWHM DiametersIntensity, ConMito, dist, Int, fitG,ContourSegment,meshSegment, distO,c, MitoIntInt_av] = measurefitFWHMROIprofileSNR(OriginalImage,pos,subNormals,mesh,meshSpacing, ContnextXsnake, ContnextYsnake,extent,Elength, StartFrame, EndFrame, segmentedMitoIDImageOld,MidPoint)
%generateFinalContours produces a final mitochondria contour and 2 contours for
%comparison based on the SIM image

%   Detailed explanation goes here
%script that generates a mitochondrial contour, calculates
%envelope curvature along contour and measures distances across mito
%according to the intensity histograms
% tic;


h = waitbar(0,'Measuring intensity profiles...');
Image=cell(1,EndFrame);
I=zeros([size(OriginalImage{StartFrame}) (EndFrame)]);
ContourSegment=cell(1,EndFrame);
meshSegment=cell(1,EndFrame);
MidPoints=cell(1,EndFrame);
for f=StartFrame:EndFrame
    Image{f} = OriginalImage{f};
  % Image{f}= im2uint16(Image{f});
     
    I(:,:,f)=Image{f};
    %Icrop(:,:,f)=imcrop(I(:,:,f),(pos./3)); 
    
    subNormals{:,f}(all(subNormals{:,f}==0,2), :)=[];  
%     TotDist{:,f}=sqrt((subNormals{:,f}(:,1)-subNormals{:,f}(:,3)).*(subNormals{:,f}(:,1)-subNormals{:,f}(:,3))+(subNormals{:,f}(:,2)-subNormals{:,f}(:,4)).*(subNormals{:,f}(:,2)-subNormals{:,f}(:,4)));
    
%     subNormals=cellfun(@(x) x*px,subNormals1,'un',0);

%     p
waitbar((f-StartFrame)/(EndFrame-StartFrame));
    if Elength<size(extent{f},1)
        if extent{f}(Elength,1)>size(subNormals{1,f},1) 
            ContourSegment{1,f}=[subNormals{1,f}(end:extent{f}(Elength,2),1) subNormals{1,f}(end:extent{f}(Elength,2),2) subNormals{1,f}(end:extent{f}(Elength,2),3) subNormals{1,f}(end:extent{f}(Elength,2),4)];
            meshSegment{1,f}=[mesh{1,f}(end:extent{f}(Elength,2),1) mesh{1,f}(end:extent{f}(Elength,2),2) mesh{1,f}(end:extent{f}(Elength,2),3) mesh{1,f}(end:extent{f}(Elength,2),4)];
            MidPoints{1,f}=[MidPoint{1,f}(end:extent{f}(Elength,2),1) MidPoint{1,f}(end:extent{f}(Elength,2),2)];
        elseif extent{f}(Elength,2)>size(subNormals{1,f},1)
            ContourSegment{1,f}=[subNormals{1,f}(extent{f}(Elength,1):end,1) subNormals{1,f}(extent{f}(Elength,1):end,2) subNormals{1,f}(extent{f}(Elength,1):end,3) subNormals{1,f}(extent{f}(Elength,1):end,4)];
            meshSegment{1,f}=[mesh{1,f}(extent{f}(Elength,1):end,1) mesh{1,f}(extent{f}(Elength,1):end,2) mesh{1,f}(extent{f}(Elength,1):end,3) mesh{1,f}(extent{f}(Elength,1):end,4)];
            MidPoints{1,f}=[MidPoint{1,f}(extent{f}(Elength,1):end,1) MidPoint{1,f}(extent{f}(Elength,1):end,2)];
        else
            ContourSegment{1,f}=[subNormals{1,f}(extent{f}(Elength,1):extent{f}(Elength,2),1) subNormals{1,f}(extent{f}(Elength,1):extent{f}(Elength,2),2) subNormals{1,f}(extent{f}(Elength,1):extent{f}(Elength,2),3) subNormals{1,f}(extent{f}(Elength,1):extent{f}(Elength,2),4)];
            meshSegment{1,f}=[mesh{1,f}(extent{f}(Elength,1):extent{f}(Elength,2),1) mesh{1,f}(extent{f}(Elength,1):extent{f}(Elength,2),2) mesh{1,f}(extent{f}(Elength,1):extent{f}(Elength,2),3) mesh{1,f}(extent{f}(Elength,1):extent{f}(Elength,2),4)];
            MidPoints{1,f}=[MidPoint{1,f}(extent{f}(Elength,1):extent{f}(Elength,2),1) MidPoint{1,f}(extent{f}(Elength,1):extent{f}(Elength,2),2)];
        end
    else
        ContourSegment{1,f}=[subNormals{1,f}(extent{f}(end,1):extent{f}(end,2),1) subNormals{1,f}(extent{f}(end,1):extent{f}(end,2),2) subNormals{1,f}(extent{f}(end,1):extent{f}(end,2),3) subNormals{1,f}(extent{f}(end,1):extent{f}(end,2),4)];
        meshSegment{1,f}=[mesh{1,f}(extent{f}(end,1):extent{f}(end,2),1) mesh{1,f}(extent{f}(end,1):extent{f}(end,2),2) mesh{1,f}(extent{f}(end,1):extent{f}(end,2),3) mesh{1,f}(extent{f}(end,1):extent{f}(end,2),4)];
        MidPoints{1,f}=[MidPoint{1,f}(extent{f}(end,1):extent{f}(end,2),1) MidPoint{1,f}(extent{f}(end,1):extent{f}(end,2),2)];
    end
    if f==StartFrame
        MaxContSegLength=length(ContourSegment{1,f});
    else
        if length(ContourSegment{1,f})>MaxContSegLength
            MaxContSegLength=length(ContourSegment{1,f});
        end
    end
    
    parfor j=1:length(ContourSegment{1,f})
        [cxtest, cytest, ctest, xitest, yitest]=improfile(I(:,:,f), [ContourSegment{1,f}(j,1), ContourSegment{1,f}(j,3)], [ContourSegment{1,f}(j,2),ContourSegment{1,f}(j,4)], 'bilinear' );
        cx{j,f}=cxtest;
        cy{j,f}=cytest;
        c{j,f}=ctest;
        xi{j,f}=xitest;
        yi{j,f}=yitest;
    end
end
close(h);
 
[SignalIdx, NoiseIdx]=measurefitFWHM_SNR(StartFrame, EndFrame, segmentedMitoIDImageOld, cx, cy, ContourSegment, MaxContSegLength);

Int=cell(MaxContSegLength,(EndFrame));
Dx=cell(MaxContSegLength,(EndFrame));
Dy=cell(MaxContSegLength,(EndFrame));
distO=cell(MaxContSegLength,(EndFrame));
dist=cell(MaxContSegLength,(EndFrame));
IntInt=cell(1,(EndFrame));
radii=cell(1,(EndFrame));
subsegmentArea=cell(1,(EndFrame));
subsegmentVol=cell(1,(EndFrame));
RatioR=cell(1,(EndFrame));
pks=cell(MaxContSegLength,(EndFrame));
locs=cell(MaxContSegLength,(EndFrame));
width=cell(MaxContSegLength,(EndFrame));
prominence=cell(MaxContSegLength,(EndFrame));
IDXpks=cell(MaxContSegLength,(EndFrame));
IDXpksS=cell(MaxContSegLength,(EndFrame));
IDXpksN=cell(MaxContSegLength,(EndFrame));
valpks=cell(MaxContSegLength,(EndFrame));
vallocs=cell(MaxContSegLength,(EndFrame));
valwidth=cell(MaxContSegLength,(EndFrame));
valprominence=cell(MaxContSegLength,(EndFrame));
IDXval=cell(MaxContSegLength,(EndFrame));
IDXvalS=cell(MaxContSegLength,(EndFrame));
IDXvalN=cell(MaxContSegLength,(EndFrame));
pkdist=cell(MaxContSegLength,(EndFrame));
Maxrow=cell(MaxContSegLength,(EndFrame));
PkInt=cell(MaxContSegLength,(EndFrame));
Pkdist=cell(MaxContSegLength,(EndFrame));
PkID=cell(MaxContSegLength,(EndFrame));
VHup=zeros(MaxContSegLength,(EndFrame));
VHlow=zeros(MaxContSegLength,(EndFrame));
VH=zeros(MaxContSegLength,(EndFrame));
PHup=zeros(MaxContSegLength,(EndFrame));
PHlow=zeros(MaxContSegLength,(EndFrame));
PH=zeros(MaxContSegLength,(EndFrame));
IDXvalR=cell(MaxContSegLength,(EndFrame));
IDXvalL=cell(MaxContSegLength,(EndFrame));
IDXfit=cell(MaxContSegLength,(EndFrame));
IDXpksR=cell(MaxContSegLength,(EndFrame));
IDXpksL=cell(MaxContSegLength,(EndFrame));
IDXfilN=cell(MaxContSegLength,(EndFrame));
fitG=cell(MaxContSegLength,(EndFrame));
Param=cell(MaxContSegLength,(EndFrame));
FWHM=cell(MaxContSegLength,(EndFrame));
CentreFit=cell(MaxContSegLength,(EndFrame));
IntO=cell(MaxContSegLength,(EndFrame));
IDX2=cell(MaxContSegLength,(EndFrame));
h = waitbar(0,'Finding valleys and peaks...');

for i=StartFrame:EndFrame
    waitbar((i-StartFrame)/(EndFrame-StartFrame));

    for j=1:length(ContourSegment{1,i})
        Int{j,i}=c{j,i};
        if min(Int{j,i})<0
%         MinI{j,i}=min(c{j,i});
            Int{j,i}=c{j,i}-min(c{j,i});
        else
            Int{j,i}=c{j,i}; %background subtracted profile
        end
        Int{j,i}(isnan(Int{j,i}(:,1)),:)=[];
        % Int{j,i}=lssmooth(Int{j,i},4);

        %transform real space coordinates to distance along line and compute
        %integrated intensity for each profile
        Dx{j,i}=cx{j,i}-xi{j,i}(1,1);
        Dy{j,i}=cy{j,i}-yi{j,i}(1,1);
        distO{j,i}= sqrt((Dx{j,i}.*Dx{j,i})+(Dy{j,i}.*Dy{j,i})); %this is the x distance along x-axis
        dist{j,i}= sqrt((Dx{j,i}.*Dx{j,i})+(Dy{j,i}.*Dy{j,i})); %this is the x distance along x-axis
        dist{j,i}(isnan(dist{j,i}(:,1)),:)=[];
        
%         if size(Int{j,i},1)==size(dist{j,i},1);
            IntInt{1,i}(j,:)=trapz(dist{j,i}, Int{j,i}); %integrated intensity for each profile
%         end

        %calculate truncated cone area
        radii{1,i}=sqrt((meshSegment{1,i}(:,1)-meshSegment{1,i}(:,3)).^2+(meshSegment{1,i}(:,2)-meshSegment{1,i}(:,4)).^2)/2;
        
        for r=1:(length(meshSegment{1,i})-1)
            subsegmentArea{1,i}(r,:)=(pi()*(radii{1,i}(r,1)+radii{1,i}(r+1,1))*sqrt((radii{1,i}(r,1)-radii{1,i}(r+1,1))^2+(meshSpacing)^2));
            subsegmentVol{1,i}(r, :)=(pi()*(1/3)*(meshSpacing)*((radii{1,i}(r,1))^2+(radii{1,i}(r+1,1))^2)+radii{1,i}(r+1,1)*radii{1,i}(r,1));
            RatioR{1,i}(r, :)=radii{1,i}(r+1,1)/radii{1,i}(r,1);
        end

        % find peaks 
        [pks{j,i},locs{j,i},width{j,i}, prominence{j,i}] = findpeaks(Int{j,i},dist{j,i}); %pks and locs are the Y and X values of the peaks
        [~,IDXpks{j,i},~]=intersect(dist{j,i},locs{j,i});

            % find peaks in SIGNAL
            IDXpksS{j,i}=intersect(IDXpks{j,i},SignalIdx{j,i});

            % find peaks in NOISE
            IDXpksN{j,i}=intersect(IDXpks{j,i},NoiseIdx{j,i});

        % find valleys
        [valpks{j,i},vallocs{j,i},valwidth{j,i}, valprominence{j,i}] = findpeaks(-Int{j,i},dist{j,i}); %pks and locs are the Y and X values of the peaks
        [~, IDXval{j,i}, ~]=intersect(dist{j,i},vallocs{j,i});
        
            % find valleys in SIGNAL
            IDXvalS{j,i}=intersect(IDXval{j,i},SignalIdx{j,i});

            % find valleys in NOISE
            IDXvalN{j,i}=intersect(IDXval{j,i},NoiseIdx{j,i});
            
       % find CENTRAL peak
       
            % if signal peaks exist
            if isempty(IDXpksS{j,i})==0
                % find main peak %MidPoint
                d1=sqrt((ContourSegment{1,i}(j,3)-MidPoints{1,i}(j,1))^2 + (ContourSegment{1,i}(j,4)-MidPoints{1,i}(j,2))^2);
                d=sqrt((ContourSegment{1,i}(j,1)-ContourSegment{1,i}(j,3))^2 + (ContourSegment{1,i}(j,2)-ContourSegment{1,i}(j,4))^2);
                pkdist{j,i}=(d1/d)*max(dist{j,i});
                knnsearch(locs{j,i},pkdist{j,i},'k',1);
                Maxrow{1,i}(:,j)=knnsearch(locs{j,i},pkdist{j,i},'k',1);
                PkInt{j,i}=pks{j,i}(Maxrow{1,i}(:,j));
                Pkdist{j,i}=locs{j,i}(Maxrow{1,i}(:,j));
%                 PkID{j,i}=knnsearch(dist{j,i},pkdist{j,i},'k',1);

                PkID{j,i}=find(dist{j,i}==locs{j,i}(Maxrow{1,i}(:,j)));

%                 PkInt{j,i}=???
%                 Pkdist{j,i}=???
%                 PkID{j,i}=???
                
            % find valley height boundaries
                if isempty(IDXvalS{j,i})==0
                    VHup(j,i)=max(Int{j,i}(IDXvalS{j,i}))/PkInt{j,i};
                else
                    VHup(j,i)=1;
                end
                if isempty(IDXvalN{j,i})==0
%                     Int{j,i}(IDXvalN{j,i})
%                     max(Int{j,i}(IDXvalN{j,i}))
%                     PkInt{j,i}
                    VHlow(j,i)=max(Int{j,i}(IDXvalN{j,i}))/PkInt{j,i};
                else
                    VHlow(j,i)=0;
                end
                if VHup(j,i)>VHlow(j,i)
                    VH(j,i)=(VHup(j,i)+VHlow(j,i))/2;
                else
                    VH(j,i)=0.99*min([VHup(j,i) VHlow(j,i)]);
                end
            
            % find peak height boundaries
                if isempty(IDXpksS{j,i})==0
                    PHlow(j,i)=max(Int{j,i}(IDXpksS{j,i}))/PkInt{j,i};
                else
                    PHlow(j,i)=1;
                end
%                 if isempty(IDXpksN{j,i})==0
%                     PHup(j,i)=max(Int{j,i}(IDXpksN{j,i}))/PkInt{j,i};
%                 else
                    PHup(j,i)=max(Int{j,i}(NoiseIdx{j,i}))/PkInt{j,i};
%                 end
                if PHup(j,i)>PHlow(j,i)
                    PH(j,i)=(PHup(j,i)+PHlow(j,i))/2;
                else
                    PH(j,i)=1.01*max([PHup(j,i) PHlow(j,i)]);
                end
            else
                PkInt{j,i}=max(Int{j,i});
                Pkdist{j,i}=mean(dist{j,i});
                PkID{j,i}=round(mean(1:length(dist{j,i})));
                PH(j,i)=1.01*max(Int{j,i})/PkInt{j,i};
                VH(j,i)=0;
            end
%             [VH(j,i) PH(j,i)]
        % cut away the noise valleys
            
        IDXvalR{j,i}=IDXvalN{j,i}(dist{j,i}(IDXvalN{j,i})<dist{j,i}(PkID{j,i}));
        IDXvalL{j,i}=IDXvalN{j,i}(dist{j,i}(IDXvalN{j,i})>dist{j,i}(PkID{j,i}));
        
        IDXvalR{j,i}=IDXvalR{j,i}(Int{j,i}(IDXvalR{j,i})<VH(j,i)*PkInt{j,i});
        IDXvalL{j,i}=IDXvalL{j,i}(Int{j,i}(IDXvalL{j,i})<VH(j,i)*PkInt{j,i});
        IDXvalR{j,i}=max(IDXvalR{j,i});
        IDXvalL{j,i}=min(IDXvalL{j,i});
        
        if isempty(IDXvalR{j,i})==1
            IDXvalR{j,i}=1;
        end
        if isempty(IDXvalR{j,i})==1
            IDXvalL{j,i}=length(dist{j,i});
        end
        
        if IDXvalR{j,i}<IDXvalL{j,i}
            IDXfit{j,i}=[IDXvalR{j,i}:IDXvalL{j,i}];
        else
            IDXfit{j,i}=SignalIdx{j,i};
        end
        
        % cut away the noise peaks
                
        IDXpksR{j,i}=IDXpksN{j,i}(dist{j,i}(IDXpksN{j,i})<dist{j,i}(PkID{j,i}));
        IDXpksL{j,i}=IDXpksN{j,i}(dist{j,i}(IDXpksN{j,i})>dist{j,i}(PkID{j,i}));            
        
        IDXpksR{j,i}=IDXpksR{j,i}(Int{j,i}(IDXpksR{j,i})>PH(j,i)*PkInt{j,i});
        IDXpksL{j,i}=IDXpksL{j,i}(Int{j,i}(IDXpksL{j,i})>PH(j,i)*PkInt{j,i});
        IDXpksR{j,i}=max(IDXpksR{j,i});
        IDXpksL{j,i}=min(IDXpksL{j,i});
        
        if isempty(IDXpksR{j,i})==1
            IDXpksR{j,i}=1;
        end
        if isempty(IDXpksL{j,i})==1
            IDXpksL{j,i}=length(dist{j,i});
        end
        
        if IDXpksR{j,i}<IDXpksL{j,i}
            if IDXpksR{j,i}>=min(IDXfit{j,i})
                IDXpksR{j,i}=min(IDXval{j,i}(IDXval{j,i}>IDXpksR{j,i}));
                IDXfit{j,i}=[IDXpksR{j,i}:max(IDXfit{j,i})];
            end
            if IDXpksL{j,i}<=max(IDXfit{j,i})
                IDXpksL{j,i}=max(IDXval{j,i}(IDXval{j,i}<IDXpksL{j,i}));
                IDXfit{j,i}=[min(IDXfit{j,i}):IDXpksL{j,i}];
            end
        end
        if isempty(IDXfit{j,i})==1
            if min(SignalIdx{j,i})>1
                minIDXfit=min(SignalIdx{j,i})-1;
            else
                minIDXfit=1;
            end
            if max(SignalIdx{j,i})<length(dist{j,i})
                maxIDXfit=max(SignalIdx{j,i})+1;
            else
                maxIDXfit=length(dist{j,i});
            end
            IDXfit{j,i}=[minIDXfit:maxIDXfit];
        end
        
        % background correction
        IDXfilN{j,i}=intersect(IDXfit{j,i},NoiseIdx{j,i});
        if isempty(NoiseIdx{j,i})==0
            Int{j,i}=Int{j,i}-median(Int{j,i}(NoiseIdx{j,i}));
        else
            Int{j,i}=Int{j,i}-min(Int{j,i}(IDXfit{j,i}));
        end
        
        dist{j,i}=dist{j,i}(IDXfit{j,i});
        Int{j,i}=Int{j,i}(IDXfit{j,i});
        
        
%     end
%     for j=1:length(ContourSegment{1,i})
        
        % fit with gaussian
%         [fitGtest Paramtest FWHMtest CentreFittest] = tryfitGaussianSNR(dist{j,i},Int{j,i},ContourSegment{1,i}(j,:),NoiseIdx{j,i}, distO{j,i},c{j,i},OriginalImage{i})
%         fitG{j,i}=fitGtest;
%         Param{j,i}=Paramtest;
%         FWHM{j,i}=FWHMtest;
%         CentreFit{j,i}=CentreFittest;

% fit with gaussian
        trywhole=1;
        Err=1;
        while Err==1;
            try
                fitG{j,i}=fit(dist{j,i}, (Int{j,i}), 'gauss1');
                Param{j,i}=coeffvalues(fitG{j,i});
                FWHM{j,i}=2.354*(1/sqrt(2))*Param{j,i}(3);
                CentreFit{j,i}=Param{j,i}(2);
%                 Total=Total+1;
                Err=0;

                if FWHM{j,i}>(1.2*(sqrt((ContourSegment{1,i}(j,1)-ContourSegment{1,i}(j,3))^2 + (ContourSegment{1,i}(j,2)-ContourSegment{1,i}(j,4))^2)))
                    ME=error('too big');
                    Err=1;
                end
                    
            catch ME
%                 display(ME);
                    if trywhole==1
                            Err=1;
                            trywhole=0;
                            dist{j,i}=distO{j,i};
                            Int{j,i}=c{j,i};
                            if isempty(NoiseIdx{j,i})==0
                                Int{j,i}=Int{j,i}-median(Int{j,i}(NoiseIdx{j,i}));
                            else
                                Int{j,i}=Int{j,i}-min(Int{j,i});
                            end
                    else
                         Err=1;
                         display(sprintf('Error in frame %d segment %d: please select profile', i,j));
                         Fig=figure;
                         subplot(1,2,1)
                         imshow(OriginalImage{i},[])
                         hold on
                         line([ContourSegment{1,i}(j,1), ContourSegment{1,i}(j,3)], [ContourSegment{1,i}(j,2),ContourSegment{1,i}(j,4)]);
                         subplot(1,2,2)
                         plot(distO{j,i},(c{j,i}-min(c{j,i})))
                         title('Select left cutoff position')
                        %      datacursormode on
                        %      dcm_obj = datacursormode;
                        [xL,yL]=ginput(1);
                        title('Select right cutoff position')
                        [xR yR]=ginput(1);
                        IntO{j,i}=(c{j,i}-min(c{j,i}));
                        Int{j,i}=IntO{j,i};
                        dist{j,i}=distO{j,i};
                        IDX2{j,i}=[1:length(dist{j,i})];
                        IDX2{j,i}(dist{j,i}<xL)=[];
                        Int{j,i}(dist{j,i}<xL)=[];
                        dist{j,i}(dist{j,i}<xL)=[];
                        Int{j,i}(dist{j,i}>xR)=[];
                        IDX2{j,i}(dist{j,i}>xR)=[];
                        dist{j,i}(dist{j,i}>xR)=[];
                        
                        if isempty(NoiseIdx{j,i})==0
                            Int{j,i}=Int{j,i}-median(IntO{j,i}(NoiseIdx{j,i}));
                        else
                            Int{j,i}=Int{j,i}-min(Int{j,i});
                        end
                        close(Fig)
                        clear ME
                    end
            end
        end
    end
end
close(h)


%% for plotting
xLeft=cell(MaxContSegLength,EndFrame);
xRight=cell(MaxContSegLength,EndFrame);
DiametersFWHM=cell(1,(EndFrame));
delta_x=cell(MaxContSegLength,EndFrame);
delta_y=cell(MaxContSegLength,EndFrame);
delta_yx=cell(MaxContSegLength,EndFrame);
delta=cell(MaxContSegLength,EndFrame);
theta=cell(MaxContSegLength,EndFrame);
alpha=cell(MaxContSegLength,EndFrame);
x_shift1=cell(MaxContSegLength,EndFrame);
y_shift1=cell(MaxContSegLength,EndFrame);
x_shift2=cell(MaxContSegLength,EndFrame);
y_shift2=cell(MaxContSegLength,EndFrame);
ConMito=cell(1,(EndFrame));
for i=StartFrame:EndFrame
for j=1:length(ContourSegment{1,i})
    
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
ContourFWHMside1=cell(1,(EndFrame));
ContourFWHMside2=cell(1,(EndFrame));
for i=StartFrame:EndFrame
 ContourFWHMside1{1,i} = lssmooth([ConMito{1,i}(:,1),ConMito{1,i}(:,2)], 10);
 ContourFWHMside2{1,i} = lssmooth([ConMito{1,i}(:,3),ConMito{1,i}(:,4)], 10);
end
 IntInt_av=cell(1,(EndFrame));
 IntDensity=cell(1,(EndFrame));
 IntDensityVol=cell(1,(EndFrame));
 MitoIntInt_av=zeros(1,(EndFrame));
 %calculate intensity density
 for i=StartFrame:EndFrame
 for r=1:(length(IntInt{1,i})-1) %%warning: matrices from these loops, might be shifted from the mesh indices
%      if isempty(nonzeros(mesh{1,i}(1)))==1
IntInt_av{1,i}(r,:)= (IntInt{1,i}(r,1)+IntInt{1,i}(r+1,1))/2;
IntDensity{1,i}(r,:)= IntInt_av{1,i}(r,:)/subsegmentArea{1,i}(r,:);
IntDensityVol{1,i}(r,:)= IntInt_av{1,i}(r,:)/subsegmentVol{1,i}(r,:);
     
 end
 
MitoIntInt_av(1,i)=sum(IntInt_av{1,i}); %sum integrated intensity for full mito
% end
 end
 
 LocInf=cell(1,(EndFrame));
ZeroLoc=cell(1,(EndFrame));
MinLoc=cell(1,(EndFrame));
MinVal=cell(1,(EndFrame));
LenRadii=cell(1,(EndFrame));

for ss=StartFrame:EndFrame
         
LocInf{ss}=find(IntDensity{1,ss}==Inf);
IntDensity{1,ss}(LocInf{ss})=[];
IntDensityVol{1,ss}(LocInf{ss})=[];

ZeroLoc{ss}=find(radii{ss}==0);%avoid detection of diameter=0
radii{ss}(ZeroLoc{ss})=[];

LenRadii{ss}=length(radii{ss});
[MinVal{ss}, MinLoc{ss}]=min(radii{ss}(5:(LenRadii{ss}-5))); %to do: make user select a boundary over which to find the min

        
end
BleachedIntDen=zeros(1,(EndFrame));
thVol=cell(1,(EndFrame));
Radius_Sol=cell(1,(EndFrame));
DiametersIntensity=cell(1,(EndFrame));
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



for rr=1:(length(IntInt_av{1,f}))
    Radius_Sol{1,f}(rr, :)=sqrt(thVol{1,f}(rr, :)*3*(1/meshSpacing)*(1/pi())*(1+(RatioR{1,f}(rr,:))^2+RatioR{1,f}(rr,:))^(-1));
    DiametersIntensity{1,f}(rr,:)=2* Radius_Sol{1,f}(rr, :);
end
end

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
        
        
 
 
%  toc;