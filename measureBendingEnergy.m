% close all; clear all;
% load('160211_comp143_event4');

function [energy, energyDensity, averageCurvatures, segmentRadii, IDmin, curvature, side1curvatures, side2curvatures, Extent, submeshCurvature, energyInner, energyDensityInner, meshspacing,segLength, smallestID, subsegment, subsegmentArea, subsegmentInner, subsegmentAreaInner]=measureBendingEnergy(StartFrame, EndFrame, ContnextXsnake, ContnextYsnake, mesh, k1, IDmin, subNormals, X_pixShift, meshSpacing, PixelSize, SaveFolder, Timescale, Delta)

%% TO DO
% solve hemispherical cap curvature calculation
% color-code mito to local energy density
bbone=cell(1,(EndFrame));
bboneLine_start=cell(1,(EndFrame));
bboneLine_finish=cell(1,(EndFrame));
X_start=cell(1,(EndFrame));
Y_start=cell(1,(EndFrame));
StartFit=cell(1,(EndFrame));
Slope_start=cell(1,(EndFrame));
YIntercept_start=cell(1,(EndFrame));
X_finish=cell(1,(EndFrame));
Y_finish=cell(1,(EndFrame));
FinishFit=cell(1,(EndFrame));
Slope_finish=cell(1,(EndFrame));
YIntercept_finish=cell(1,(EndFrame));
Start=cell(1,(EndFrame));
Finish=cell(1,(EndFrame));
NewXpoint_start=cell(1,(EndFrame));
NewYpoint_start=cell(1,(EndFrame));
NewXpoint_finish=cell(1,(EndFrame));
NewYpoint_finish=cell(1,(EndFrame));
NewPoint_start=cell(1,(EndFrame));
NewPoint_finish=cell(1,(EndFrame));
polygonBorder=cell(1,(EndFrame));
in=cell(1,(EndFrame));
Side1=cell(1,(EndFrame));
ks1=cell(1,(EndFrame));
Side2=cell(1,(EndFrame));
ks2=cell(1,(EndFrame));
curvature=cell(2,(EndFrame));
%generate the backbone matrix
for f=StartFrame:EndFrame

for l=1:length(subNormals{1,f})
      bbone{f}(l,:)=[0.5*(subNormals{1,f}(l,1)+subNormals{1,f}(l,3)) 0.5*(subNormals{1,f}(l,2)+subNormals{1,f}(l,4))];
end
Len=length(subNormals{1,f});
bboneLine_start{f}=bbone{f}(1:4,:); % first 4 points of backbone line for frame f
bboneLine_finish{f}=bbone{f}(Len-4:Len, :); %last 4 points of backbone line for frame f
end

%extend the backbone line in each direction
h = waitbar(0,'Calculating bending energy...');

for t=StartFrame:EndFrame
                waitbar((t-StartFrame)/(EndFrame-StartFrame));

%start of the backbone
   X_start{t}=bboneLine_start{t}(1, 1)-bboneLine_start{t}(4, 1); %displacement in x direction of top/start of backbone
   Y_start{t}=bboneLine_start{t}(1, 2)-bboneLine_start{t}(4, 2); %displacement in y direction of top/start of backbone
   StartFit{t}=polyfit(bboneLine_start{t}(:,1),bboneLine_start{t}(:,2), 1);
   Slope_start{t}=StartFit{t}(1,1);
   YIntercept_start{t}=StartFit{t}(1,2);
   
 %end of the backbone
   X_finish{t}=bboneLine_finish{t}(4, 1)-bboneLine_finish{t}(1, 1); %displacement in x direction of top/start of backbone
   Y_finish{t}=bboneLine_finish{t}(4, 2)-bboneLine_finish{t}(1, 2); %displacement in y direction of top/start of backbone
   FinishFit{t}=polyfit(bboneLine_finish{t}(:,1),bboneLine_finish{t}(:,2), 1);
   Slope_finish{t}=FinishFit{t}(1,1);
   YIntercept_finish{t}=FinishFit{t}(1,2);
   
   
   
   Start{t}=[X_start{t} Y_start{t}];
   
   Finish{t}=[X_finish{t} Y_finish{t}];
   
   %define new point from start of bbone line
   
   if sign(Start{t})==[1,1]; %case 1, x,y positive displacement
       
       NewXpoint_start{t}=bboneLine_start{t}(1, 1)+X_pixShift;
       NewYpoint_start{t}=NewXpoint_start{t}*Slope_start{t}+YIntercept_start{t};
   end
       
   if sign(Start{t})==[1,0]
       NewXpoint_start{t}=bboneLine_start{t}(1, 1)+X_pixShift;
       NewYpoint_start{t}=NewXpoint_start{t}*Slope_start{t}+YIntercept_start{t};
   end
   
   if sign(Start{t})==[0,1]
       NewXpoint_start{t}=bboneLine_start{t}(1, 1)+X_pixShift;
       NewYpoint_start{t}=NewXpoint_start{t}*Slope_start{t}+YIntercept_start{t};
   end
   if sign(Start{t})==[1,-1];%case 2, x positive, y negative displacement
       
       NewXpoint_start{t}=bboneLine_start{t}(1, 1)+X_pixShift;
       NewYpoint_start{t}=NewXpoint_start{t}*Slope_start{t}+YIntercept_start{t};
   end
   
   
    if sign(Start{t})==[-1,-1]; %case 3, x,y negative displacement
       
       NewXpoint_start{t}=bboneLine_start{t}(1, 1)-X_pixShift;
       NewYpoint_start{t}=NewXpoint_start{t}*Slope_start{t}+YIntercept_start{t};
    end
   
    if sign(Start{t})==[-1,1]; %case 4, x negative, y positive displacement
       
       NewXpoint_start{t}=bboneLine_start{t}(1, 1)-X_pixShift;
       NewYpoint_start{t}=NewXpoint_start{t}*Slope_start{t}+YIntercept_start{t};
    end
    
    if sign(Start{t})==[-1,0]
       NewXpoint_start{t}=bboneLine_start{t}(1, 1)-X_pixShift;
       NewYpoint_start{t}=NewXpoint_start{t}*Slope_start{t}+YIntercept_start{t};
    end
    
    if sign(Start{t})==[0,-1]
       NewXpoint_start{t}=bboneLine_start{t}(1, 1)-X_pixShift;
       NewYpoint_start{t}=NewXpoint_start{t}*Slope_start{t}+YIntercept_start{t};
    end
    
    %define new point from end of bbone line
    
    if sign(Finish{t})==[1,1]; %case 1, x,y positive displacement
       
       NewXpoint_finish{t}=bboneLine_finish{t}(4, 1)+X_pixShift;
       NewYpoint_finish{t}=NewXpoint_finish{t}*Slope_finish{t}+YIntercept_finish{t};
   end
       
   
   if sign(Finish{t})==[1,-1];%case 2, x positive, y negative displacement
       
       NewXpoint_finish{t}=bboneLine_finish{t}(4, 1)+X_pixShift;
       NewYpoint_finish{t}=NewXpoint_finish{t}*Slope_finish{t}+YIntercept_finish{t};
   end
   
   
    if sign(Finish{t})==[-1,-1]; %case 3, x,y negative displacement
       
       NewXpoint_finish{t}=bboneLine_finish{t}(4, 1)-X_pixShift;
       NewYpoint_finish{t}=NewXpoint_finish{t}*Slope_finish{t}+YIntercept_finish{t};
    end
   
    if sign(Finish{t})==[-1,1]; %case 4, x negative, y positive displacement
       
       NewXpoint_finish{t}=bboneLine_finish{t}(4, 1)-X_pixShift;
       NewYpoint_finish{t}=NewXpoint_finish{t}*Slope_finish{t}+YIntercept_finish{t};
    end
%     t
    NewPoint_start{t}=[NewXpoint_start{t} NewYpoint_start{t}];
    NewPoint_finish{t}=[NewXpoint_finish{t} NewYpoint_finish{t}];

    polygonBorder{t}=[NewPoint_start{t}; bbone{t}; NewPoint_finish{t};flipud(subNormals{t}(:, 1:2));NewPoint_start{t}];
    in{t}=inpolygon(ContnextXsnake{t}, ContnextYsnake{t}, polygonBorder{t}(:,1), polygonBorder{t}(:,2));
    
    for m=1:length(in{t})
        if in{t}(m)==1
            Side1{t}(m, :)=[ContnextXsnake{t}(m, 1) ContnextYsnake{t}(m, 1)];
            ks1{1,t}(m, :)=k1{t}(m,:);
        else
            Side2{t}(m, :)=[ContnextXsnake{t}(m, 1) ContnextYsnake{t}(m, 1)];
            ks2{1,t}(m,:)=k1{t}(m,:);
        end
    end
 Side1{t}(~any(Side1{t},2), :)=[];
 Side2{t}(~any(Side2{t},2), :)=[];
 ks1{t}(~any(ks1{t},2), :)=[];
 ks2{t}(~any(ks2{t},2), :)=[];
  
   curvature{1,t}=[Side1{t} ks1{1,t}];
   curvature{2,t}=[Side2{t} ks2{1,t}];
   
%    meshs{t}=union(meshPoints{1,t},meshPoints{2,t});

end

close(h);
%curvature is a matrix with col1: xid (mesh1,2), col2: yid(mesh3, 4), col3:
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

% separate contour and envelope curvature into 2 sides
submesh=cell(1,(EndFrame));
radii=cell(2,(EndFrame));
submeshCurvature=cell(2,(EndFrame));
segLength=cell(2,(EndFrame));
% curvside1=cell((EndFrame),);
% curvside2=cell((EndFrame));
averageRadii=cell(1,(EndFrame));
segmentRadii=cell(1,(EndFrame));
smallestID=zeros(1,(EndFrame));
meshspacing=cell(1,(EndFrame));
averageCurvatures=cell(1,(EndFrame));
largestID1=zeros(1,(EndFrame));
largestID2=zeros(1,(EndFrame));
Extent=cell(1,(EndFrame));
% StartFrame=1;
for f=StartFrame:EndFrame;%size(mesh,2)
    for s=1:2 %sides 1 and 2
        submesh{f}=mesh{1,f};
        radii{s,f}=sqrt((submesh{1,f}(:,1)-submesh{1,f}(:,3)).^2+(submesh{1,f}(:,2)-submesh{1,f}(:,4)).^2)*PixelSize/2;
       
        %find closest curvature point for each mesh point
        IDX=knnsearch(curvature{s,f}(:,1:2),submesh{1,f}(:,(2*s-1):(2*s))); 
       
        for j=1:(length(IDX)-1)
            if IDX(j)>IDX(j+1)
                submeshCurvature{s,f}(j,1)=mean(curvature{s,f}(IDX(j+1):IDX(j),3));
                submeshCurvature{s,f}(j,2)=1/((radii{s,f}(j)+radii{s,f}(j+1))/2);
                segLength{s,f}(j)=sum(sqrt(sum((diff(curvature{s,f}(IDX(j+1):IDX(j),1:2)).^2),2)));
                if s==1
                    curvside1{f,j}=curvature{s,f}(IDX(j+1):IDX(j),3);
                elseif s==2
                    curvside2{f,j}=curvature{s,f}(IDX(j+1):IDX(j),3);
                end

            elseif IDX(j)<IDX(j+1)
                submeshCurvature{s,f}(j,1)=mean(curvature{s,f}(IDX(j):IDX(j+1),3));
                submeshCurvature{s,f}(j,2)=1/((radii{s,f}(j)+radii{s,f}(j+1))/2);
                segLength{s,f}(j)=sum(sqrt(sum((diff(curvature{s,f}(IDX(j):IDX(j+1),1:2)).^2),2)));

                if s==1
                    curvside1{f,j}=curvature{s,f}(IDX(j):IDX(j+1),3);
                elseif s==2
                    curvside2{f,j}=curvature{s,f}(IDX(j):IDX(j+1),3);
                end
            else
                submeshCurvature{s,f}(j,1)=0;
                submeshCurvature{s,f}(j,2)=0;
                segLength{s,f}(j)=0;
                if s==1
                    curvside1{f,j}=0;
                elseif s==2
                    curvside2{f,j}=0;
                end
               
            end
        end
    end
    
        averageRadii{f}=(radii{1,f}+radii{2,f})/2;
        segmentRadii{f}=(averageRadii{1,f}(1:(end-1))+averageRadii{1,f}(2:end))/2;
        smallestID(f)=IDmin(f);%find(segmentRadii{1,f}(meshs{1,f},:)==min(segmentRadii{1,f}(meshs{1,f})));
%         IDXmesh=knnsearch(submesh{1,f},submesh{1,f}(meshs{1,f}(smallestID(f)),:));
        smallestID(f)=IDmin(f);%IDXmesh;
        while smallestID(f)>size(segmentRadii{f},1)
            smallestID(f)=IDmin(f)-1;
        end
            
        if f==StartFrame
            absSmallestID=f;
        else
            if segmentRadii{1,f}(smallestID(f))<segmentRadii{1,absSmallestID}(smallestID(absSmallestID))
                absSmallestID=f;
            end
        end
        
    for i=1:(size(submesh{f},1)-1)
%         [bbone{f}(i,:); bbone{f}(i+1,:)]
        meshspacing{f}(i)=sum(sqrt(sum((diff([bbone{f}(i,:);bbone{f}(i+1,:)]).^2),2)));
%         sum(sqrt(sum((diff([ContnextXsnake{f}(IDX(j):IDX(j+1),1) ContnextYsnake{f}(IDX(j):IDX(j+1),1)]).^2),2)));
% [bbone{f}(i+1,:);bbone{f}(i)]
        if segLength{1,f}(i)>(10*meshspacing{f}(i))
            segLength{1,f}(i)=10*meshspacing{f}(i);
        end
        if segLength{2,f}(i)>(10*meshspacing{f}(i))
            segLength{2,f}(i)=10*meshspacing{f}(i);
        end
        averageCurvatures{f}(i,1)=(submeshCurvature{1,f}(i,1)*(segLength{1,f}(i)/(segLength{1,f}(i)+segLength{2,f}(i)))+submeshCurvature{2,f}(i,1)*((segLength{2,f}(i)/(segLength{1,f}(i)+segLength{2,f}(i)))));
        averageCurvatures{f}(i,2)=(submeshCurvature{1,f}(i,2)+submeshCurvature{2,f}(i,2))/2;
        subsegment(f,i)=(averageCurvatures{1,f}(i,1)+averageCurvatures{1,f}(i,2))^2 * ...
            (pi()*(averageRadii{1,f}(i)+averageRadii{1,f}(i+1))*sqrt((averageRadii{1,f}(i)-averageRadii{1,f}(i+1))^2+(meshspacing{f}(i))^2));
        subsegmentArea(f,i)=(pi()*(averageRadii{1,f}(i)+averageRadii{1,f}(i+1))*sqrt((averageRadii{1,f}(i)-averageRadii{1,f}(i+1))^2+(meshspacing{f}(i)*PixelSize)^2));
        subsegmentInner(f,i)=(((averageCurvatures{1,f}(i,1)).^(-1)-Delta).^(-1)+((averageCurvatures{1,f}(i,2)).^(-1)-Delta).^(-1))^2 * ...
            (pi()*(averageRadii{1,f}(i)-Delta+averageRadii{1,f}(i+1)-Delta)*sqrt((averageRadii{1,f}(i)-Delta-averageRadii{1,f}(i+1)+Delta)^2+(meshspacing{f}(i)*PixelSize)^2));
        subsegmentAreaInner(f,i)=(pi()*(averageRadii{1,f}(i)-Delta+averageRadii{1,f}(i+1)-Delta)*sqrt((averageRadii{1,f}(i)-Delta-averageRadii{1,f}(i+1)+Delta)^2+(meshspacing{f}(i)*PixelSize)^2));
    end
    
    % find local maxima
    largestID2(f)=find(segmentRadii{1,f}((smallestID(f):end),:)==max(segmentRadii{1,f}((smallestID(f):end),:)))+smallestID(f)-1;
    largestID1(f)=find(segmentRadii{1,f}((1:smallestID(f)),:)==max(segmentRadii{1,f}((1:smallestID(f)),:)));

    for i=(smallestID(f)):(size(segmentRadii{1,f},1)-1)
        if MaxMaxBandPc==1
            if(segmentRadii{1,f}(i)>=segmentRadii{1,f}(largestID2(f),:)*BandPc)
                localMaxID2(f)=i;
                break;
            end
        elseif MaxMaxBandFix==1
            if(segmentRadii{1,f}(i)>=segmentRadii{1,f}(largestID2(f),:)-BandFix)
                localMaxID2(f)=i;
                break;
            end
        else
            if(segmentRadii{1,f}(i)>segmentRadii{1,f}(i+1))
                
                localMaxID2(f)=i;
                break;
            end
        end
    end
    for i=(smallestID(f)):-1:2
        if MaxMaxBandPc==1
            if(segmentRadii{1,f}(i)>=segmentRadii{1,f}(largestID1(f),:)*BandPc)
                localMaxID1(f)=i;
                break;
            end
        elseif MaxMaxBandFix==1
            if(segmentRadii{1,f}(i)>=segmentRadii{1,f}(largestID1(f),:)-BandFix)
                localMaxID1(f)=i;
                break;
            end
        else
            if(segmentRadii{1,f}(i)>segmentRadii{1,f}(i-1))
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
                energy(f,1)=sum(subsegment(f,localMaxID1(f):localMaxID2(f)));
                area(f,1)=sum(subsegmentArea(f,localMaxID1(f):localMaxID2(f)));
                energyDensity(f,1)=energy(f,1)/area(f,1);
                Length(f,1)=sum(sqrt(diff(bbone{1,f}(localMaxID1(f):localMaxID2(f),1)).^2 + diff(bbone{1,f}(localMaxID1(f):localMaxID2(f),2)).^2));
                frame(f,1)=f;
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
                energy(f,seg+1)=sum(subsegment(f,lowerLimit:upperLimit));
                area(f,seg+1)=sum(subsegmentArea(f,lowerLimit:upperLimit));
                energyDensity(f,seg+1)=energy(f,seg+1)/area(f,seg+1);
                Length(f,seg+1)=sum(sqrt(diff(bbone{1,f}(lowerLimit:upperLimit,1)).^2 + diff(bbone{1,f}(lowerLimit:upperLimit,2)).^2));
                frame(f,seg+1)=f;
                Extent{f}(seg+1,1)=lowerLimit;
                Extent{f}(seg+1,2)=upperLimit;
                
                energyInner(f,seg+1)=sum(subsegmentInner(f,lowerLimit:upperLimit));
                areaInner(f,seg+1)=sum(subsegmentAreaInner(f,lowerLimit:upperLimit));
                energyDensityInner(f,seg+1)=energyInner(f,seg+1)/areaInner(f,seg+1);
                                seg=seg+1;


        end
    end
    
    
    
    
%     energy(f)=((submeshCurvature{1,f}(:,1)+submeshCurvature{2,f}(:,1))/2 + (submeshCurvature{1,f}(:,2)+submeshCurvature{2,f}(:,2))/2).^2 * ...
%         (pi()*((radii{1,f}(1:(end-1))+radii{2,f}(1:(end-1)))/2 + (radii{1,f}(2:end)+radii{2,f}(2:end))/2)*...
%         sqrt(((radii{1,f}(1:(end-1))+radii{2,f}(1:(end-1)))/2 - (radii{1,f}(2:end)+radii{2,f}(2:end))/2).^2+(meshSpacing*PixelSize)^2));
%     energy(f)=sum((((submeshCurvature{1,f}(:,2)+submeshCurvature{2,f}(:,2))./2+(submeshCurvature{1,f}(:,1)+submeshCurvature{2,f}(:,1))./2).^2)*2*pi()*(sqrt((submesh{1,f}(j,1)-submesh{1,f}(j,3)).^2+(submesh{1,f}(j,2)-submesh{1,f}(j,4))^2)*30/2*1));
%     energy(f)=sum((((submeshCurvature{1,f}(:,2)+submeshCurvature{2,f}(:,2))./2+(submeshCurvature{1,f}(:,1)+submeshCurvature{2,f}(:,1))./2).^2)*2*pi()*(sqrt((submesh{1,f}(j,1)-submesh{1,f}(j,3)).^2+(submesh{1,f}(j,2)-submesh{1,f}(j,4))^2)*30/2*1))/...
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
        
            energy(f,1)=sum(subsegment(f,localMaxID1(absSmallestID):localMaxID2(absSmallestID)));
            area(f,1)=sum(subsegmentArea(f,localMaxID1(absSmallestID):localMaxID2(absSmallestID)));
            energyDensity(f,1)=energy(f,1)/area(f,1);
            Length(f,1)=sum(sqrt(diff(bbone{1,f}(localMaxID1(absSmallestID):localMaxID2(absSmallestID),1)).^2 + diff(bbone{1,f}(localMaxID1(absSmallestID):localMaxID2(absSmallestID),2)).^2));
            frame(f,1)=f;
        
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

% Zen=[reshape(Length,[prod(size(Length)) 1]) reshape(frame,[prod(size(frame)) 1]) reshape(energy,[prod(size(energy)) 1])];
% ZenD=[reshape(Length,[prod(size(Length)) 1]) reshape(frame,[prod(size(frame)) 1]) reshape(energyDensity,[prod(size(energyDensity)) 1])];


if MaxMax==1
%     figure
%     plot((-size(energy,1):-1)*Timescale,energy(:,1)')
%     title('Energy');
%     xlabel('Frame number')
%     ylabel('Elastic energy [a.u.]')
%    savefig([SaveFolder '\BEplot']);
% 
%     figure
%     plot((-size(energyDensity,1):-1)*Timescale,energyDensity(:,1)')
%     title('Energy Density');
%     xlabel('Frame number')
%     ylabel('Elastic energy density [a.u.]')
%    savefig([SaveFolder '\BEDplot']);

else
    
%     figure
%     plot((-size(energy,1):-1)*Timescale,energy(:,segSizes)')
%     title('Energy outer membrane');
%     xlabel('Frame number')
%     ylabel('Bending energy [a.u.]')
%    savefig([SaveFolder '\BEplot']);
% 
%     figure
%     plot((-size(energyDensity,1):-1)*Timescale,energyDensity(:,segSizes)')
%     title('Energy Density outer membrane');
%     xlabel('Frame number')
%     ylabel('Bending energy density [a.u.]')
%    savefig([SaveFolder '\BEDplot']);
%     
%     figure
%     plot((-size(energyInner,1):-1)*Timescale,energyInner(:,segSizes+1)')
%     title('Energy inner membrane');
%     xlabel('Frame number')
%     ylabel('Bending energy [a.u.]')
%    savefig([SaveFolder '\BEplotI']);
% 
%     figure
%     plot((-size(energyDensityInner,1):-1)*Timescale,energyDensityInner(:,segSizes+1)')
%     title('Energy Density inner membrane');
%     xlabel('Frame number')
%     ylabel('Bending energy density [a.u.]')
%    savefig([SaveFolder '\BEDplotI']);

end

if MaxMax==1
    
else
%     figure
%     surf(energy)
%     ylabel('Frame number')
%     xlabel('Segment length [pixel]')
%     zlabel('Bending energy [a.u.]')
%     colormap hot
%    savefig([SaveFolder '\BEsurf']);
% 
%     figure
%     surf(energyDensity)
%     ylabel('Frame number')
%     xlabel('Segment length [pixel]')
%     zlabel('Bending energy density [a.u.]')
%     colormap hot
%    savefig([SaveFolder '\BEDsurf']);

end
