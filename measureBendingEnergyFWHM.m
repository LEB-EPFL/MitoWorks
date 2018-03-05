% close all; clear all;
% load('160211_comp143_event4');

function [energyFWHM, energyDensityFWHM, averageCurvaturesFWHM, segmentRadiiFWHM, IDX_fwhm, curvatureFWHM, Extent, energyFWHMInner, energyDensityFWHMInner, side1curvatures, side2curvatures, meshspacing, segLength, smallestID, subsegment, subsegmentArea, subsegmentInner, subsegmentAreaInner]=measureBendingEnergyFWHM(StartFrame, EndFrame, ContourFWHMside1, ContourFWHMside2, ks1,ks2, IDX_fwhm,PixelSize, SaveFolder,Delta)
%% TO DO
% solve hemispherical cap curvatureFWHM calculation
% color-code mito to local energyFWHM density


%generate the backbone matrix
% for f=StartFrame:EndFrame
% 
% for l=1:length(subNormals{1,f})
%       bbone{f}(l,:)=[0.5*(subNormals{1,f}(l,1)+subNormals{1,f}(l,3)) 0.5*(subNormals{1,f}(l,2)+subNormals{1,f}(l,4))];
% end
% Len=length(subNormals{1,f});
% bboneLine_start{f}=bbone{f}(1:4,:); % first 4 points of backbone line for frame f
% bboneLine_finish{f}=bbone{f}(Len-4:Len, :); %last 4 points of backbone line for frame f
% end
% 
% %extend the backbone line in each direction
% h = waitbar(0,'Calculating bending energyFWHM...');
h = waitbar(0,'Calculating bending energy for FWHM...');
bbone=cell(1,(EndFrame));
% bboneLine_start=cell(1,(EndFrame));
% bboneLine_finish=cell(1,(EndFrame));
% X_start=cell(1,(EndFrame));
% Y_start=cell(1,(EndFrame));
% StartFit=cell(1,(EndFrame));
% Slope_start=cell(1,(EndFrame));
% YIntercept_start=cell(1,(EndFrame));
% X_finish=cell(1,(EndFrame));
% Y_finish=cell(1,(EndFrame));
% FinishFit=cell(1,(EndFrame));
% Slope_finish=cell(1,(EndFrame));
% YIntercept_finish=cell(1,(EndFrame));
% Start=cell(1,(EndFrame));
% Finish=cell(1,(EndFrame));
% NewXpoint_start=cell(1,(EndFrame));
% NewYpoint_start=cell(1,(EndFrame));
% NewXpoint_finish=cell(1,(EndFrame));
% NewYpoint_finish=cell(1,(EndFrame));
% NewPoint_start=cell(1,(EndFrame));
% NewPoint_finish=cell(1,(EndFrame));
% polygonBorder=cell(1,(EndFrame));
% in=cell(1,(EndFrame));
% Side1=cell(1,(EndFrame));
% ks1=cell(1,(EndFrame));
% Side2=cell(1,(EndFrame));
% ks2=cell(1,(EndFrame));
curvature=cell(2,(EndFrame));
curvatureFWHM=cell(2,(EndFrame));
for t=StartFrame:EndFrame
                waitbar((t-StartFrame)/(EndFrame-StartFrame));
% 
% %start of the backbone
%    X_start{t}=bboneLine_start{t}(1, 1)-bboneLine_start{t}(4, 1); %displacement in x direction of top/start of backbone
%    Y_start{t}=bboneLine_start{t}(1, 2)-bboneLine_start{t}(4, 2); %displacement in y direction of top/start of backbone
%    StartFit{t}=polyfit(bboneLine_start{t}(:,1),bboneLine_start{t}(:,2), 1);
%    Slope_start{t}=StartFit{t}(1,1);
%    YIntercept_start{t}=StartFit{t}(1,2);
%    
%  %end of the backbone
%    X_finish{t}=bboneLine_finish{t}(4, 1)-bboneLine_finish{t}(1, 1); %displacement in x direction of top/start of backbone
%    Y_finish{t}=bboneLine_finish{t}(4, 2)-bboneLine_finish{t}(1, 2); %displacement in y direction of top/start of backbone
%    FinishFit{t}=polyfit(bboneLine_finish{t}(:,1),bboneLine_finish{t}(:,2), 1);
%    Slope_finish{t}=FinishFit{t}(1,1);
%    YIntercept_finish{t}=FinishFit{t}(1,2);
%    
%    
%    
%    Start{t}=[X_start{t} Y_start{t}];
%    
%    Finish{t}=[X_finish{t} Y_finish{t}];
%    
%    %define new point from start of bbone line
%    
%    if sign(Start{t})==[1,1]; %case 1, x,y positive displacement
%        
%        NewXpoint_start{t}=bboneLine_start{t}(1, 1)+X_pixShift;
%        NewYpoint_start{t}=NewXpoint_start{t}*Slope_start{t}+YIntercept_start{t};
%    end
%        
%    
%    if sign(Start{t})==[1,-1];%case 2, x positive, y negative displacement
%        
%        NewXpoint_start{t}=bboneLine_start{t}(1, 1)+X_pixShift;
%        NewYpoint_start{t}=NewXpoint_start{t}*Slope_start{t}+YIntercept_start{t};
%    end
%    
%    
%     if sign(Start{t})==[-1,-1]; %case 3, x,y negative displacement
%        
%        NewXpoint_start{t}=bboneLine_start{t}(1, 1)-X_pixShift;
%        NewYpoint_start{t}=NewXpoint_start{t}*Slope_start{t}+YIntercept_start{t};
%     end
%    
%     if sign(Start{t})==[-1,1]; %case 4, x negative, y positive displacement
%        
%        NewXpoint_start{t}=bboneLine_start{t}(1, 1)-X_pixShift;
%        NewYpoint_start{t}=NewXpoint_start{t}*Slope_start{t}+YIntercept_start{t};
%     end
%     
%     %define new point from end of bbone line
%     
%     if sign(Finish{t})==[1,1]; %case 1, x,y positive displacement
%        
%        NewXpoint_finish{t}=bboneLine_finish{t}(4, 1)+X_pixShift;
%        NewYpoint_finish{t}=NewXpoint_finish{t}*Slope_finish{t}+YIntercept_finish{t};
%    end
%        
%    
%    if sign(Finish{t})==[1,-1];%case 2, x positive, y negative displacement
%        
%        NewXpoint_finish{t}=bboneLine_finish{t}(4, 1)+X_pixShift;
%        NewYpoint_finish{t}=NewXpoint_finish{t}*Slope_finish{t}+YIntercept_finish{t};
%    end
%    
%    
%     if sign(Finish{t})==[-1,-1]; %case 3, x,y negative displacement
%        
%        NewXpoint_finish{t}=bboneLine_finish{t}(4, 1)-X_pixShift;
%        NewYpoint_finish{t}=NewXpoint_finish{t}*Slope_finish{t}+YIntercept_finish{t};
%     end
%    
%     if sign(Finish{t})==[-1,1]; %case 4, x negative, y positive displacement
%        
%        NewXpoint_finish{t}=bboneLine_finish{t}(4, 1)-X_pixShift;
%        NewYpoint_finish{t}=NewXpoint_finish{t}*Slope_finish{t}+YIntercept_finish{t};
%     end
%     NewPoint_start{t}=[NewXpoint_start{t} NewYpoint_start{t}];
%     NewPoint_finish{t}=[NewXpoint_finish{t} NewYpoint_finish{t}];
% 
%     polygonBorder{t}=[NewPoint_start{t}; bbone{t}; NewPoint_finish{t};flipud(subNormals{t}(:, 1:2));NewPoint_start{t}];
%     in{t}=inpolygon(ContnextXsnake{t}, ContnextYsnake{t}, polygonBorder{t}(:,1), polygonBorder{t}(:,2));
    
%     for m=1:length(in{t})
%         if in{t}(m)==1
%             Side1{t}(m, :)=[ContnextXsnake{t}(m, 1) ContnextYsnake{t}(m, 1)];
%             ks1{1,t}(m, :)=k1{t}(m,:);
%         else
%             Side2{t}(m, :)=[ContnextXsnake{t}(m, 1) ContnextYsnake{t}(m, 1)];
%             ks2{1,t}(m,:)=k1{t}(m,:);
%         end
%     end
%  Side1{t}(~any(Side1{t},2), :)=[];
%  Side2{t}(~any(Side2{t},2), :)=[];
%  ks1{t}(~any(ks1{t},2), :)=[];
%  ks2{t}(~any(ks2{t},2), :)=[];
  
   curvatureFWHM{1,t}=[ContourFWHMside1{t} ks1{1,t}];
   curvatureFWHM{2,t}=[ContourFWHMside2{t} ks2{1,t}];
   
   curvature{1,t}=[ContourFWHMside1{t} ks1{1,t}];
   curvature{2,t}=[ContourFWHMside2{t} ks2{1,t}];
   
%    meshs{t}=union(meshPoints{1,t},meshPoints{2,t});

end

close(h);
%curvatureFWHM is a matrix with col1: xid (mesh1,2), col2: yid(mesh3, 4), col3:
%k1 (to be split into2)
WindowSize=7;
% Offset=20/3;
% meshSpacing=1;

segSizes=[1 2 3 4 5];
segSize=10;
MaxMax=0;
MaxMaxBandPc=0; BandPc=0.9; 
MaxMaxBandFix=0; BandFix=20; 

lastFrame=0;

% MaxMax for each frame = MaxMax=1, MaxMaxBand=0, lastFrame=0
% MaxMaxBand for each frame = MaxMax=1, MaxMaxBand=1, lastFrame=0
% MaxMax from last frame = MaxMax=1, MaxMaxBand=0, lastFrame=1
% MaxMaxBand from last frame = MaxMax=1, MaxMaxBand=1, lastFrame=1

% separate contour and envelope curvatureFWHM into 2 sides
submesh=cell(1,(EndFrame));
radii=cell(2,(EndFrame));
submeshCurvature=cell(2,(EndFrame));
segLength=cell(2,(EndFrame));
% curvside1=cell((EndFrame),);
% curvside2=cell((EndFrame));
averageRadii=cell(1,(EndFrame));
averageRadiiIn=cell(1,(EndFrame));
% segmentRadii=cell(1,(EndFrame));
smallestID=zeros(1,(EndFrame));
meshspacing=cell(1,(EndFrame));
% averageCurvatures=cell(1,(EndFrame));
largestID1=zeros(1,(EndFrame));
largestID2=zeros(1,(EndFrame));
mesh=cell(1,(EndFrame));
segmentRadiiFWHM=cell(1,(EndFrame));
segmentRadiiFWHMIn=cell(1,(EndFrame));
averageCurvaturesFWHM=cell(1,(EndFrame));
averageCurvaturesFWHMIn=cell(1,(EndFrame));
Extent=cell(1,(EndFrame));
% StartFrame=1;
submeshCurvatureIn=cell(2,(EndFrame));
localMaxID2=zeros(1,(EndFrame));
for f=StartFrame:EndFrame;%size(mesh,2)
    mesh{1,f}=[ContourFWHMside1{1,f} ContourFWHMside2{1,f}];
    for s=1:2 %sides 1 and 2
        submesh{f}=mesh{1,f};
        radii{s,f}=sqrt((submesh{1,f}(:,1)-submesh{1,f}(:,3)).^2+(submesh{1,f}(:,2)-submesh{1,f}(:,4)).^2)*PixelSize/2;
       
        %find closest curvatureFWHM point for each mesh point
        IDX=knnsearch(curvatureFWHM{s,f}(:,1:2),submesh{1,f}(:,(2*s-1):(2*s))); 
       
        for j=1:(length(IDX)-1)
            if IDX(j)>IDX(j+1)
                submeshCurvature{s,f}(j,1)=mean(curvatureFWHM{s,f}(IDX(j+1):IDX(j),3));
                submeshCurvature{s,f}(j,2)=1/((radii{s,f}(j)+radii{s,f}(j+1))/2);
                segLength{s,f}(j)=sum(sqrt(sum((diff(curvatureFWHM{s,f}(IDX(j+1):IDX(j),1:2)).^2),2)));
                submeshCurvatureIn{s,f}(j,1)=mean(((curvatureFWHM{s,f}(IDX(j+1):IDX(j),3).^(-1))-Delta).^(-1));
                submeshCurvatureIn{s,f}(j,2)=1/((radii{s,f}(j)-Delta+radii{s,f}(j+1)-Delta)/2);
                if s==1
                    curvside1{f,j}=curvatureFWHM{s,f}(IDX(j+1):IDX(j),3);
                elseif s==2
                    curvside2{f,j}=curvatureFWHM{s,f}(IDX(j+1):IDX(j),3);
                end

            elseif IDX(j)<IDX(j+1)
                submeshCurvature{s,f}(j,1)=mean(curvatureFWHM{s,f}(IDX(j):IDX(j+1),3));
                submeshCurvature{s,f}(j,2)=1/((radii{s,f}(j)+radii{s,f}(j+1))/2);
                segLength{s,f}(j)=sum(sqrt(sum((diff(curvatureFWHM{s,f}(IDX(j):IDX(j+1),1:2)).^2),2)));

                if s==1
                    curvside1{f,j}=curvature{s,f}(IDX(j):IDX(j+1),3);
                elseif s==2
                    curvside2{f,j}=curvature{s,f}(IDX(j):IDX(j+1),3);
                end
                submeshCurvatureIn{s,f}(j,1)=mean(((curvatureFWHM{s,f}(IDX(j):IDX(j+1),3).^(-1))-Delta).^(-1));
                submeshCurvatureIn{s,f}(j,2)=1/((radii{s,f}(j)-Delta+radii{s,f}(j+1)-Delta)/2);
            else
                submeshCurvature{s,f}(j,1)=0;
                submeshCurvature{s,f}(j,2)=0;
                segLength{s,f}(j)=0;

                if s==1
                    curvside1{f,j}=0;
                elseif s==2
                    curvside2{f,j}=0;
                end
                
                submeshCurvatureIn{s,f}(j,1)=0;
                submeshCurvatureIn{s,f}(j,2)=0;
            end
        end
    end
    
        averageRadii{f}=(radii{1,f}+radii{2,f})/2;
        averageRadiiIn{f}=(radii{1,f}-Delta+radii{2,f}-Delta)/2;
        segmentRadiiFWHM{f}=(averageRadii{1,f}(1:(end-1))+averageRadii{1,f}(2:end))/2;
        segmentRadiiFWHMIn{f}=(averageRadiiIn{1,f}(1:(end-1))+averageRadiiIn{1,f}(2:end))/2;

        smallestID(f)=IDX_fwhm(f);%find(segmentRadiiFWHM{1,f}(meshs{1,f},:)==min(segmentRadiiFWHM{1,f}(meshs{1,f})));
%         IDXmesh=knnsearch(submesh{1,f},submesh{1,f}(meshs{1,f}(smallestID(f)),:));
        smallestID(f)=IDX_fwhm(f);%IDXmesh;
        while smallestID(f)>size(segmentRadiiFWHM{f},1)
            smallestID(f)=IDX_fwhm(f)-1;
        end
            
        if f==StartFrame
            absSmallestID=f;
        else
            if segmentRadiiFWHM{1,f}(smallestID(f))<segmentRadiiFWHM{1,absSmallestID}(smallestID(absSmallestID))
                absSmallestID=f;
            end
        end
        bbone{1,f}=[(submesh{1,f}(:,1)+submesh{1,f}(:, 3))./2 (submesh{1,f}(:,2)+submesh{1,f}(:, 4))./2 ]; % find backbone

        
    for i=1:(size(submesh{f},1)-1)
        meshspacing{f}(i)=sum(sqrt(sum((diff([bbone{f}(i,:);bbone{f}(i+1,:)]).^2),2)));
if segLength{1,f}(i)>(10*meshspacing{f}(i))
            segLength{1,f}(i)=10*meshspacing{f}(i);
        end
        if segLength{2,f}(i)>(10*meshspacing{f}(i))
            segLength{2,f}(i)=10*meshspacing{f}(i);
        end
        averageCurvaturesFWHM{f}(i,1)=(submeshCurvature{1,f}(i,1)*(segLength{1,f}(i)/(segLength{1,f}(i)+segLength{2,f}(i)))+submeshCurvature{2,f}(i,1)*(segLength{2,f}(i)/(segLength{1,f}(i)+segLength{2,f}(i))));
        averageCurvaturesFWHM{f}(i,2)=(submeshCurvature{1,f}(i,2)+submeshCurvature{2,f}(i,2))/2;
        averageCurvaturesFWHMIn{f}(i,1)=(submeshCurvatureIn{1,f}(i,1)+submeshCurvatureIn{2,f}(i,1))/2;
        averageCurvaturesFWHMIn{f}(i,2)=(submeshCurvatureIn{1,f}(i,2)+submeshCurvatureIn{2,f}(i,2))/2;
        subsegment(f,i)=(averageCurvaturesFWHM{1,f}(i,1)+averageCurvaturesFWHM{1,f}(i,2))^2 * ...
            (pi()*(averageRadii{1,f}(i)+averageRadii{1,f}(i+1))*sqrt((averageRadii{1,f}(i)-averageRadii{1,f}(i+1))^2+(meshspacing{f}(i)*PixelSize)^2));
        subsegmentArea(f,i)=(pi()*(averageRadii{1,f}(i)+averageRadii{1,f}(i+1))*sqrt((averageRadii{1,f}(i)-averageRadii{1,f}(i+1))^2+(meshspacing{f}(i)*PixelSize)^2));
        subsegmentInner(f,i)=(averageCurvaturesFWHMIn{1,f}(i,1)+averageCurvaturesFWHMIn{1,f}(i,2))^2 * ...
            (pi()*(averageRadiiIn{1,f}(i)+averageRadiiIn{1,f}(i+1))*sqrt((averageRadiiIn{1,f}(i)-averageRadiiIn{1,f}(i+1))^2+(meshspacing{f}(i)*PixelSize)^2));
        subsegmentAreaInner(f,i)=(pi()*(averageRadiiIn{1,f}(i)+averageRadiiIn{1,f}(i+1))*sqrt((averageRadiiIn{1,f}(i)-averageRadiiIn{1,f}(i+1))^2+(meshspacing{f}(i)*PixelSize)^2));
        
    end
    
    % find local maxima
    largestID2(f)=find(segmentRadiiFWHM{1,f}((smallestID(f):end),:)==max(segmentRadiiFWHM{1,f}((smallestID(f):end),:)))+smallestID(f)-1;
    largestID1(f)=find(segmentRadiiFWHM{1,f}((1:smallestID(f)),:)==max(segmentRadiiFWHM{1,f}((1:smallestID(f)),:)));

    for i=(smallestID(f)):(size(segmentRadiiFWHM{1,f},1)-1)
        if MaxMaxBandPc==1
            if(segmentRadiiFWHM{1,f}(i)>=segmentRadiiFWHM{1,f}(largestID2(f),:)*BandPc)
                localMaxID2(f)=i;
                break;
            end
        elseif MaxMaxBandFix==1
            if(segmentRadiiFWHM{1,f}(i)>=segmentRadiiFWHM{1,f}(largestID2(f),:)-BandFix)
                localMaxID2(f)=i;
                break;
            end
        else

            if(segmentRadiiFWHM{1,f}(i)>segmentRadiiFWHM{1,f}(i+1))
                
                localMaxID2(f)=i;
                break;
            end
        end
    end
    for i=(smallestID(f)):-1:2
        if MaxMaxBandPc==1
            if(segmentRadiiFWHM{1,f}(i)>=segmentRadiiFWHM{1,f}(largestID1(f),:)*BandPc)
                localMaxID1(f)=i;
                break;
            end
        elseif MaxMaxBandFix==1
            if(segmentRadiiFWHM{1,f}(i)>=segmentRadiiFWHM{1,f}(largestID1(f),:)-BandFix)
                localMaxID1(f)=i;
                break;
            end
        else
            if(segmentRadiiFWHM{1,f}(i)>segmentRadiiFWHM{1,f}(i-1))
                localMaxID1(f)=i;
                break;
            end
        end
    end 
    
    
    seg=0;
    upLimMet=0;
    lowLimMet=0;
                    bbone{1,f}=[(submesh{1,f}(:,1)+submesh{1,f}(:, 3))./2 (submesh{1,f}(:,2)+submesh{1,f}(:, 4))./2 ]; % find backbone

    if MaxMax==1
        if lastFrame==0
                energyFWHM(f,1)=sum(subsegment(f,localMaxID1(f):localMaxID2(f)));
                area(f,1)=sum(subsegmentArea(f,localMaxID1(f):localMaxID2(f)));
                energyDensityFWHM(f,1)=energyFWHM(f,1)/area(f,1);
                Length(f,1)=sum(sqrt(diff(bbone{1,f}(localMaxID1(f):localMaxID2(f),1)).^2 + diff(bbone{1,f}(localMaxID1(f):localMaxID2(f),2)).^2));
                frame(f,1)=f;
                
                energyFWHMInner(f,1)=sum(subsegmentInner(f,localMaxID1(f):localMaxID2(f)));
                areaInner(f,1)=sum(subsegmentAreaInner(f,localMaxID1(f):localMaxID2(f)));
                energyDensityFWHMInner(f,1)=energyFWHMInner(f,1)/areaInner(feg+1);
        end
    else
        while(upLimMet==0 || lowLimMet==0)
            if (seg+smallestID(f))>=(size(submesh{f},1)-1)
                upperLimit=(size(submesh{f},1)-1);
                upLimMet=1;
            else 
                upperLimit=seg+smallestID(f);
            end

            if (smallestID(f)-seg)<=1
                lowerLimit=1;
                lowLimMet=1;
            else
                lowerLimit=smallestID(f)-seg;
            end
                bbone{1,f}=[(submesh{1,f}(:,1)+submesh{1,f}(:, 3))./2 (submesh{1,f}(:,2)+submesh{1,f}(:, 4))./2 ]; % find backbone

                energyFWHM(f,seg+1)=sum(subsegment(f,lowerLimit:upperLimit));
                Extent{f}(seg+1,:)=[lowerLimit upperLimit];
                area(f,seg+1)=sum(subsegmentArea(f,lowerLimit:upperLimit));
                energyDensityFWHM(f,seg+1)=energyFWHM(f,seg+1)/area(f,seg+1);
                Length(f,seg+1)=sum(sqrt(diff(bbone{1,f}(lowerLimit:upperLimit,1)).^2 + diff(bbone{1,f}(lowerLimit:upperLimit,2)).^2));
                frame(f,seg+1)=f;
                
                energyFWHMInner(f,seg+1)=sum(subsegmentInner(f,lowerLimit:upperLimit));
                areaInner(f,seg+1)=sum(subsegmentAreaInner(f,lowerLimit:upperLimit));
                energyDensityFWHMInner(f,seg+1)=energyFWHMInner(f,seg+1)/areaInner(f,seg+1);
                
                for i=lowerLimit:upperLimit
                    if i==lowerLimit
                        side1curvatures{f,seg+1}=[]; 
                    end
                    side1curvatures{f,seg+1}=[side1curvatures{f,seg+1} curvside1{f,i}'];
                end
                for i=lowerLimit:upperLimit
                    if i==lowerLimit
                        side2curvatures{f,seg+1}=[]; 
                    end
                    side2curvatures{f,seg+1}=[side2curvatures{f,seg+1} curvside2{f,i}'];
                end
                seg=seg+1;


        end
    end
    
    
    
    
%     energyFWHM(f)=((submeshCurvature{1,f}(:,1)+submeshCurvature{2,f}(:,1))/2 + (submeshCurvature{1,f}(:,2)+submeshCurvature{2,f}(:,2))/2).^2 * ...
%         (pi()*((radii{1,f}(1:(end-1))+radii{2,f}(1:(end-1)))/2 + (radii{1,f}(2:end)+radii{2,f}(2:end))/2)*...
%         sqrt(((radii{1,f}(1:(end-1))+radii{2,f}(1:(end-1)))/2 - (radii{1,f}(2:end)+radii{2,f}(2:end))/2).^2+(meshSpacing*PixelSize)^2));
%     energyFWHM(f)=sum((((submeshCurvature{1,f}(:,2)+submeshCurvature{2,f}(:,2))./2+(submeshCurvature{1,f}(:,1)+submeshCurvature{2,f}(:,1))./2).^2)*2*pi()*(sqrt((submesh{1,f}(j,1)-submesh{1,f}(j,3)).^2+(submesh{1,f}(j,2)-submesh{1,f}(j,4))^2)*30/2*1));
%     energyFWHM(f)=sum((((submeshCurvature{1,f}(:,2)+submeshCurvature{2,f}(:,2))./2+(submeshCurvature{1,f}(:,1)+submeshCurvature{2,f}(:,1))./2).^2)*2*pi()*(sqrt((submesh{1,f}(j,1)-submesh{1,f}(j,3)).^2+(submesh{1,f}(j,2)-submesh{1,f}(j,4))^2)*30/2*1))/...
%         sum(2*pi()*(sqrt((submesh{1,f}(j,1)-submesh{1,f}(j,3)).^2+(submesh{1,f}(j,2)-submesh{1,f}(j,4))^2)*30/2*1));
    if f==StartFrame
        minMeshSize=size(mesh{1,f},1);
    else
        if size(mesh{1,f},1)<minMeshSize
            minMeshSize=size(mesh{1,f},1);
        end
    end
    
    if MaxMax==1
        if lastFrame==0
            figure
                hold on;

            plot(bbone{1,f}(localMaxID1(f):localMaxID2(f),1), bbone{1,f}(localMaxID1(f):localMaxID2(f),2), 'k-');
            plot([mesh{1,f}(:,1:2:3)],[mesh{1,f}(:,2:2:4)])
            title(sprintf('Frame %d',f));
            hold on
            for i=localMaxID1(f):localMaxID2(f)
            line([mesh{1,f}(i,1:2:3)],[mesh{1,f}(i,2:2:4)])
            end
            axis equal
            hold off
        end
    else
%         figure
%             hold on;
%         if (smallestID(f)-segSize)<=0
%     
%             plot(bbone{1,f}(1:(smallestID(f)+segSize),1), bbone{1,f}(1:(smallestID(f)+segSize),2), 'k-');
%             plot([mesh{1,f}(:,1:2:3)],[mesh{1,f}(:,2:2:4)])
%             title(sprintf('Frame %d',f));
%             hold on
%             for i=1:(smallestID(f)+segSize)
%             line([mesh{1,f}(i,1:2:3)],[mesh{1,f}(i,2:2:4)])
%             end
%             axis equal
%             hold off
%             
%         elseif (smallestID(f)+segSize)>size(bbone{1,f},1)
%             plot(bbone{1,f}((smallestID(f)-segSize):end,1), bbone{1,f}((smallestID(f)-segSize):end,2), 'k-');
%             plot([mesh{1,f}(:,1:2:3)],[mesh{1,f}(:,2:2:4)])
%             title(sprintf('Frame %d',f));
%             hold on
%             for i=(smallestID(f)-segSize):size(mesh{1,f},1)
%             line([mesh{1,f}(i,1:2:3)],[mesh{1,f}(i,2:2:4)])
%             end
%             axis equal
%             hold off
%         else
%             plot(bbone{1,f}((smallestID(f)-segSize):(smallestID(f)+segSize),1), bbone{1,f}((smallestID(f)-segSize):(smallestID(f)+segSize),2), 'k-');
%             plot([mesh{1,f}(:,1:2:3)],[mesh{1,f}(:,2:2:4)])
%             title(sprintf('Frame %d',f));
%             hold on
%             for i=(smallestID(f)-segSize):(smallestID(f)+segSize)
%             line([mesh{1,f}(i,1:2:3)],[mesh{1,f}(i,2:2:4)])
%             end
%             axis equal
%             hold off
%         end
    end
end

if MaxMax==1 && lastFrame==1
%     absSmallestID
    for f=1:size(mesh,2)
        
            energyFWHM(f,1)=sum(subsegment(f,localMaxID1(absSmallestID):localMaxID2(absSmallestID)));
            area(f,1)=sum(subsegmentArea(f,localMaxID1(absSmallestID):localMaxID2(absSmallestID)));
            energyDensityFWHM(f,1)=energyFWHM(f,1)/area(f,1);
            Length(f,1)=sum(sqrt(diff(bbone{1,f}(localMaxID1(absSmallestID):localMaxID2(absSmallestID),1)).^2 + diff(bbone{1,f}(localMaxID1(absSmallestID):localMaxID2(absSmallestID),2)).^2));
            frame(f,1)=f;
            
            energyFWHMInner(f,1)=sum(subsegmentInner(f,localMaxID1(absSmallestID):localMaxID2(absSmallestID)));
            areaInner(f,1)=sum(subsegmentAreaInner(f,localMaxID1(absSmallestID):localMaxID2(absSmallestID)));
            energyDensityFWHMInner(f,1)=energyFWHMInner(f,1)/areaInner(f,1);
        
            figure
            hold on;

            plot(bbone{1,f}(localMaxID1(absSmallestID):localMaxID2(absSmallestID),1), bbone{1,f}(localMaxID1(absSmallestID):localMaxID2(absSmallestID),2), 'k-');
            plot([mesh{1,f}(:,1:2:3)],[mesh{1,f}(:,2:2:4)])
            title(sprintf('Frame %d',f));
            hold on
            for i=localMaxID1(absSmallestID):localMaxID2(absSmallestID)
                line([mesh{1,f}(i,1:2:3)],[mesh{1,f}(i,2:2:4)])
            end
            axis equal
            hold off
    end
end

% Lengths=mean(Length);

% Zen=[reshape(Length,[prod(size(Length)) 1]) reshape(frame,[prod(size(frame)) 1]) reshape(energyFWHM,[prod(size(energyFWHM)) 1])];
% ZenD=[reshape(Length,[prod(size(Length)) 1]) reshape(frame,[prod(size(frame)) 1]) reshape(energyDensityFWHM,[prod(size(energyDensityFWHM)) 1])];


if MaxMax==1
 
else
    
%     figure
%     plot(1:size(energyFWHM,1),energyFWHM(:,segSizes)')
%     title('Energy outer membrane');
%     xlabel('Frame number')
%     ylabel('Bending energyFWHM [a.u.]')
%     savefig([SaveFolder '\fwhmBEplot']);

%     figure
%     plot(1:size(energyDensityFWHM,1),energyDensityFWHM(:,segSizes)')
%     title('Energy Density outer membrane');
%     xlabel('Frame number')
%     ylabel('Bending energyFWHM density [a.u.]')
%     savefig([SaveFolder '\fwhmBEDplot']);
%     
%     figure
%     plot(1:size(energyFWHMInner,1),energyFWHMInner(:,segSizes)')
%     title('Energy inner membrane');
%     xlabel('Frame number')
%     ylabel('Bending energyFWHM [a.u.]')
%     savefig([SaveFolder '\fwhmBEplotI']);
% 
%     figure
%     plot(1:size(energyDensityFWHMInner,1),energyDensityFWHMInner(:,segSizes)')
%     title('Energy Density inner membrane');
%     xlabel('Frame number')
%     ylabel('Bending energyFWHM density [a.u.]')
%     savefig([SaveFolder '\fwhmBEDplotI']);

end

if MaxMax==1
    
else
%     figure
%     surf(energyFWHM)
%     ylabel('Frame number')
%     xlabel('Segment length [pixel]')
%     zlabel('Bending energyFWHM [a.u.]')
%     colormap hot
%     savefig([SaveFolder '\fwhmBEsurf']);
% 
%     figure
%     surf(energyDensityFWHM)
%     ylabel('Frame number')
%     xlabel('Segment length [pixel]')
%     zlabel('Bending energyFWHM density [a.u.]')
%     colormap hot
%     savefig([SaveFolder '\fwhmBEDsurf']);

end


