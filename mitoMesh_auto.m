function [mesh, Contour, subNormals, DiametersSnake,MidPoint]=mitoMesh_auto(fname,ContnextX,ContnextY,PixelSize, Offset,startFrame,endFrame, pos,meshD, SpurC)
% script that creates an ordered mesh of the mitochondrion from the
% segmented image and a given contour 
%
% Inputs:
% fname = path to .tif file with the segmented image
% contourFilename = name of .mat file with the contours (ContnextX,ContnextY)
% PixelSize = pixel size [nm]
% Offset = image offset [nm]
% startFrame = start frame
% endFrame = end frame
% Tau = smoothing parameter, larger tau means smoother curve, try 25
%
% Outputs:
% mesh = mesh{f}(i,:) coordinates of mesh points frame f, line i in format
%         [x1 y1 x2 y2]
% Contour = [ContnextX ContnextY] format
% subNormals = subNormals{f}(i,:) coordinates of mesh points frame f, line
%               i in format [x1 y1 x2 y2]

% clear all; close all;
set(0,'DefaultFigureWindowStyle','docked')
%% ===== INPUT AND OUTPUT FILE NAMES =====

addpath('C:\Users\laboleb\Documents\MATLAB\lssmooth');

% input file
% fname = '160211Composite143_event4_8pt1s.tif';
image_source = 'SIM';
% output file name
Resize=1;
regroupmesh=1;
ExtensionDistance=round(200/PixelSize);

Tau=askTau();

%% START

% [info, numberOfFrames,width,height] = getImageInfo(fname);

% startFrame=1;
% endFrame=numberOfFrames;
% imshow(fname{1, startFrame});
% 

ComposedImage=(zeros(size(im2bw(fname{startFrame}))));
for f=startFrame:endFrame
    ComposedImage=max(ComposedImage,im2bw(fname{f}));
end
imshow(imresize(imfill(ComposedImage, 'holes'),Resize))
title('Drag a box around the mito.')
box=imrect;
Positions=box.getPosition();
close;
% h = waitbar(0,'Running mitoMesh.m');
FirstBranching=1;
Contour=cell(1,(endFrame));
imageBW=cell(1,(endFrame));
resized=cell(1,(endFrame));
% skelD=cell(1,(endFrame));
% ySkel=cell(1,(endFrame));
% xSkel=cell(1,(endFrame));
% ySkelNew=cell(1,(endFrame));
% xSkelNew=cell(1,(endFrame));
% xSkelNewFinal=cell(1,(endFrame));
% ySkelNewFinal=cell(1,(endFrame));
SmoothBackbone=cell(1,(endFrame));
Grad=cell(1,(endFrame));
InvGrad=cell(1,(endFrame));
MidPoint=cell(1,(endFrame));
Intercept=cell(1,(endFrame));
InvIntercept=cell(1,(endFrame));
Normals=cell(1,(endFrame));
mesh=cell(1,(endFrame));
DiametersSnake=cell(1,(endFrame));
subNormals=cell(1,(endFrame));
for f = startFrame:endFrame
	
%     waitbar((f-startFrame)/(endFrame-startFrame));
	[height, width]=size(fname{f});
    fname{f}=im2bw(fname{f});
    % try width - Contnext if the contour is flipped
   % Contour{f}=[(ContnextX{f}./PixelSize-Offset) ContnextY{f}./PixelSize-Offset];
   Contour{f}=[ContnextX{f} ContnextY{f}];
    clear B; clear E; clear yS; clear xS; clear yB; clear xB; clear yE; clear xE; clear distE; clear maxEr; clear maxEc;
%     while(Repeat=='1')
    imageBW{f} = fname{f};
    resized{f} = imresize(imfill(imageBW{f}, 'holes'),Resize);
    hold off
    skel= bwmorph(resized{f},'thin',Inf);
    skel=bwmorph(skel,'spur',SpurC);
    B = bwmorph(skel, 'branchpoints');
    E = bwmorph(skel, 'endpoints');
    
    imshow(resized{f},[]);
    axis([Positions(1) Positions(1)+Positions(3) Positions(2) Positions(2)+Positions(4)]);
    [~, MSGID] = lastwarn();
    warning('off', MSGID);
    hold all;
    [yS,xS] = find(skel); plot(xS,yS, 'c.')
    [yB,xB] = find(B);  plot(xB,yB,'gx'); %posB= [xB yB];
    [yE,xE] = find(E); plot(xE,yE,'bx'); %posE= [xE yE]; 
%     distE=pdist([yE,xE]);
%     distE=squareform(distE);
%     [maxEr maxEc]=find(distE==max(max(distE)),1);
    
    categoriesE=zeros(length(xE),1);
    
    % check edges
%     Ef=zeros(size(E));
    Extremes=zeros(size(E));
    if f>startFrame
        IDXE=[];
        OldEs=[oldxE' oldyE'];
    end
            oldxE=[]; oldyE=[];

    if length(xE)==2
        for i=1:length(xE)
            B(yE(i),xE(i))=1;
            Extremes(yE(i),xE(i))=1;
            oldyE=[oldyE yE(i)]; oldxE=[oldxE xE(i)];
        end
    else
        if f==startFrame
            oldyE=[]; oldxE=[];
        for i=1:length(xE)
            OK=0;
            while OK==0
                imshow(resized{f});
                axis([Positions(1) Positions(1)+Positions(3) Positions(2) Positions(2)+Positions(4)]);

                title('Choose end points');
                hold on
                plot(xS,yS, 'c.')
                plot(xB,yB,'gx')
                plot(xE,yE,'bx')
                plot(xE(i),yE(i),'ro')
                choiceE=input(sprintf('Keep (1) or discard (2) for frame %d?',f));
                if choiceE==1
                    B(yE(i),xE(i))=1;
                    Extremes(yE(i),xE(i))=1;
                    oldyE=[oldyE yE(i)]; oldxE=[oldxE xE(i)];
                    OK=1;
                elseif choiceE==2
                    OK=1;
                elseif choiceE==3
                    assignE=input('Assign endpoint to which subgroup (insert number)?');
                    categoriesE(i)=assignE;
                    OK=1;
                else
                end
            end
        end
        else
%             size(OldEs)
%             size([xE yE])
            IDX=knnsearch([xE yE],OldEs);
            numel(IDX)
            B(yE(IDX),xE(IDX))=1;
            Extremes(yE(IDX),xE(IDX))=1;
            oldyE=[yE(IDX)']; oldxE=[xE(IDX)'];

        end
            
    end
    


            
    
    
    % check branches
    for i=1:length(xB)
       
    end

    Dmask = false(size(skel));
    B_loc = find(B);
    for k = 1:numel(xE)
        D = bwdistgeodesic(skel,xE(k),yE(k));
        distanceToBranchPt = min(D(B_loc));
        Dmask(D < distanceToBranchPt) =true; clear distanceToBranchPt;
    end
    skelD{f} = skel - Dmask;
    hold off
    imshow(skelD{f});
    axis([Positions(1) Positions(1)+Positions(3) Positions(2) Positions(2)+Positions(4)]);

    hold all;
    [y,x] = find(B);
    plot(x,y,'ro')
    
    skelD{f}=skelD{f};
        [Ystart, Xstart]=find(Extremes);

    
    [ySkel{f},xSkel{f}]=find(skelD{f}); 
    Err=1;
    while Err==1;
    %put back to original size
    ySkel{f}=ySkel{f}/Resize; xSkel{f}=xSkel{f}/Resize; 


%   find start
        IDX=zeros(1,length(2:size(ySkel{f},1)));

    IDX(1)=knnsearch([ySkel{f} xSkel{f}],  [Ystart(1)/Resize Xstart(1)/Resize]);
    
    X=[ySkel{f} xSkel{f}];
    for i=2:size(ySkel{f},1)
        if i==2
            IDXtemp=knnsearch(X,X(IDX(i-1),:),'K', numel(IDX)+1);
            Temp=setdiff(IDXtemp,IDX);
            IDX(i)=Temp(knnsearch(X(Temp,:),X(IDX(i-1),:)));
        else
            IDXtemp=knnsearch(X,X(IDX(i-1),:),'K', numel(IDX)+1);
            Temp=setdiff(IDXtemp,IDX);
            IDX(i)=Temp(knnsearch(X(Temp,:),X(IDX(i-1),:)));
        end
    end
    clear Temp; clear IDXtemp;
    ySkel{f}=ySkel{f}(IDX);
    xSkel{f}=xSkel{f}(IDX);
    
    ySkelold=ySkel{f};
    xSkelold=xSkel{f};
    
    % pick fewer values
    Resize2=Resize*regroupmesh;
    for i=1:size(ySkel{f},1)/Resize2
        ySkelNew{f}(i)=mean(ySkelold(((i-1)*Resize2+1):(i*Resize2)));
        xSkelNew{f}(i)=mean(xSkelold(((i-1)*Resize2+1):(i*Resize2)));
    end
    
%     ySkelNewFinal{f}(1:2:(2*length(ySkelNew{f})-1))=ySkelNew{f};
%     ySkelNewFinal{f}(2:2:(2*length(ySkelNew{f})-2))=(ySkelNew{f}(1:(end-1))+ySkelNew{f}(2:end))./2;
%     xSkelNewFinal{f}(1:2:(2*length(xSkelNew{f})-1))=xSkelNew{f};
%     xSkelNewFinal{f}(2:2:(2*length(xSkelNew{f})-2))=(xSkelNew{f}(1:(end-1))+xSkelNew{f}(2:end))./2;
    
    for j=1:meshD
%         [size(j:meshD:(meshD*(length(xSkelNew{f})-1))) size(xSkelNew{f}(1:(end-1))) size(xSkelNew{f}(2:end))]
        xSkelNewFinal{f}(j:meshD:(meshD*(length(xSkelNew{f})-1)))=(((meshD-(j-1))/meshD).*xSkelNew{f}(1:(end-1)))+((((j-1))/meshD).*xSkelNew{f}(2:end));
        ySkelNewFinal{f}(j:meshD:(meshD*(length(ySkelNew{f})-1)))=(((meshD-(j-1))/meshD).*ySkelNew{f}(1:(end-1)))+((((j-1))/meshD).*ySkelNew{f}(2:end));
    end
    Err=0;
    try
    hold off
    imshow(imageBW{f});
    axis([Positions(1) Positions(1)+Positions(3) Positions(2) Positions(2)+Positions(4)]);

    hold on
    line(xSkelNewFinal{f},ySkelNewFinal{f}, 'Color', 'red');
    plot(xSkelNewFinal{f},ySkelNewFinal{f},'rx');
    
    %smooth backbone
    SmoothBackbone{f}=lssmooth([xSkelNewFinal{f}' ySkelNewFinal{f}'],str2double(Tau{1})*meshD);
    plot(SmoothBackbone{f}(:,1),SmoothBackbone{f}(:,2),'b');
    plot(Contour{f}(:,1),Contour{f}(:,2),'g', 'linewidth', 1.5)
    if isempty(SmoothBackbone{f})==1 || numel(xSkelNewFinal{f})<=3
                Err=1;
                error('err message');
    end
    catch
        Err=1;
            warning('Too few points: please drag backbone line manually.')
            BBLine=imline;
            BBBurnedImage=BBLine.createMask();
            [BBrs BBcs]=find(BBBurnedImage);
%             figure 
%             imshow(BBBurnedImage);
%             pause;
%             close;
            
            xSkel{f}=BBcs;
            ySkel{f}=BBrs;
            
            Ystart=BBrs(1);
            Xstart=BBcs(1);
            title('Please draw backbone line');
    end
    end
    NoSegments=size(SmoothBackbone{f},1)-1;
    
Intersections=cell((endFrame),length(NoSegments));
subIntersections1=cell((endFrame),length(NoSegments));
subIntersections2=cell((endFrame),length(NoSegments));

    
    for i=1:NoSegments
        Grad{f}(i)=(SmoothBackbone{f}(i+1,2)-SmoothBackbone{f}(i,2))/(SmoothBackbone{f}(i+1,1)-SmoothBackbone{f}(i,1));
        InvGrad{f}(i)=-1/(Grad{f}(i));
        MidPoint{f}(i,:)=[(SmoothBackbone{f}(i+1,1)+SmoothBackbone{f}(i,1))/2 (SmoothBackbone{f}(i+1,2)+SmoothBackbone{f}(i,2))/2];
        Intercept{f}(i)=MidPoint{f}(i,2)-Grad{f}(i)*MidPoint{f}(i,1);
%         InvIntercept{f}(i)=(Grad{f}(i)-Intercept{f}(i))*MidPoint{f}(i,1)+Intercept{f}(i);
        InvIntercept{f}(i)=MidPoint{f}(i,2)-InvGrad{f}(i)*MidPoint{f}(i,1);
%         hold on
%         refline(InvGrad{f}(i),InvIntercept{f}(i));
%         hold off
        Normals{f}(i,:)=[0 width (InvGrad{f}(i)*0+InvIntercept{f}(i)) (InvGrad{f}(i)*width+InvIntercept{f}(i))];
        Intersections{f,i}=InterX([Normals{f}(i,1) Normals{f}(i,2); Normals{f}(i,3) Normals{f}(i,4)],[Contour{f}(:,1)'; Contour{f}(:,2)']);
        subIntersections1{f,i}=Intersections{f,i}(:,Intersections{f,i}(2,:)<(Grad{f}(i)*Intersections{f,i}(1,:)+Intercept{f}(i)));
        subIntersections2{f,i}=Intersections{f,i}(:,Intersections{f,i}(2,:)>(Grad{f}(i)*Intersections{f,i}(1,:)+Intercept{f}(i)));

%     [MidPoint{f}(i,:)]
%     Intersections{f,i}
%     [subIntersections1{f,i}]
%     [subIntersections2{f,i}]
        IDX1=knnsearch([subIntersections1{f,i}]',[MidPoint{f}(i,:)],'K',1);
        IDX2=knnsearch([subIntersections2{f,i}]',[MidPoint{f}(i,:)],'K',1);
%         IDX=knnsearch([Intersections{f,i}]',[MidPoint{f}(i,:)],'K',2);
%         Intersections{f,i}=Intersections{f,i}(:,IDX);
        Intersections{f,i}=[subIntersections1{f,i}(:,IDX1) subIntersections2{f,i}(:,IDX2)];
        hold on
        if isempty(Intersections{f,i})==0
            if size(Intersections{f,i})==[2,1]
                DeltaX=sqrt((ExtensionDistance^2)/(1+InvGrad{f}(i)^2));
                DeltaY=DeltaX*InvGrad{f}(i);
                mesh{f}(i,:)=[Intersections{f,i}(1,1) Intersections{f,i}(2,1) Intersections{f,i}(1,1) Intersections{f,i}(2,1)];
                
                DiametersSnake{f}(i,:)=sqrt((Intersections{f,i}(1,1)-Intersections{f,i}(1,1))^2 + (Intersections{f,i}(2,1)-Intersections{f,i}(2,1))^2);
                subNormals{f}(i,:)=[Intersections{f,i}(1,1)+DeltaX Intersections{f,i}(2,1)+DeltaY Intersections{f,i}(1,1)-DeltaX Intersections{f,i}(2,1)-DeltaY];
            else
                line([Intersections{f,i}(1,1) Intersections{f,i}(1,2)],[Intersections{f,i}(2,1) Intersections{f,i}(2,2)]);
                mesh{f}(i,:)=[Intersections{f,i}(1,1) Intersections{f,i}(2,1) Intersections{f,i}(1,2) Intersections{f,i}(2,2)];
                DiametersSnake{f}(i,:)=sqrt((Intersections{f,i}(1,1)-Intersections{f,i}(1,2))^2 + (Intersections{f,i}(2,1)-Intersections{f,i}(2,2))^2);
                
                
                DeltaX=sqrt((ExtensionDistance^2)/(1+InvGrad{f}(i)^2));
                DeltaY=DeltaX*InvGrad{f}(i);
                
                    subNormals{f}(i,:)=[Intersections{f,i}(1,1)+DeltaX Intersections{f,i}(2,1)+DeltaY Intersections{f,i}(1,2)-DeltaX Intersections{f,i}(2,2)-DeltaY];
             
                if pdist([subNormals{f}(i,1) subNormals{f}(i,2); subNormals{f}(i,3) subNormals{f}(i,4)])<pdist([mesh{f}(i,1) mesh{f}(i,2); mesh{f}(i,3) mesh{f}(i,4)])+ExtensionDistance*1.9
                    subNormals{f}(i,:)=[Intersections{f,i}(1,1)-DeltaX Intersections{f,i}(2,1)-DeltaY Intersections{f,i}(1,2)+DeltaX Intersections{f,i}(2,2)+DeltaY];
                end
                for j=1:4
                    if subNormals{f}(i,j)<1
                        subNormals{f}(i,j)=1;
                    elseif (j==1 || j==3) && subNormals{f}(i,j)>width
                        subNormals{f}(i,j)=width;
                    elseif (j==2 || j==4) && subNormals{f}(i,j)>height
                        subNormals{f}(i,j)=height;
                    end
                end
                plot([subNormals{f}(i,1) subNormals{f}(i,3)],[subNormals{f}(i,2) subNormals{f}(i,4)],'xc');
            end
        end

        % order mesh
    end
    for i=1:length(mesh{f}(:,1))
        if i>1
            IDX=knnsearch([mesh{f}(i,1:2); mesh{f}(i,3:4)],mesh{f}(i-1,1:2));
            if IDX==2
                mesh{f}(i,:)=[mesh{f}(i,3) mesh{f}(i,4) mesh{f}(i,1) mesh{f}(i,2)] ;     
            end
        end
        
        if i>1
            IDX=knnsearch([subNormals{f}(i,1:2); subNormals{f}(i,3:4)],subNormals{f}(i-1,1:2));
            if IDX==2
                subNormals{f}(i,:)=[subNormals{f}(i,3) subNormals{f}(i,4) subNormals{f}(i,1) subNormals{f}(i,2)] ;     
            end
        end
        
        
    end
                            drawnow;

%     
%     %find if any loops
%     if (xSkelNew{f}(1)<xSkelNew{f}(2) && (isempty(find(diff(xSkelNew{f})<0))==0))
%             flipIDX(f)=min(find(diff(xSkelNew{f})<0));
%     elseif(xSkelNew{f}(1)>xSkelNew{f}(2) && isempty(find(diff(xSkelNew{f})==0)))
%             flipIDX(f)=min(find(diff(xSkelNew{f})>0));
%     else
%         flipIDX(f)=1;
%     end
%     plot(xSkelNew{f}(flipIDX(f)),ySkelNew{f}(flipIDX(f)),'bx');
%     
%     for i=1:flipIDX(f)-1
%         xvals1(1+(i-1)*xspacing:i*xspacing)=linspace(xSkelNew{f}(i),xSkelNew{f}(i+1),xspacing);
%     end
%     for i=flipIDX(f):size(xSkelNew{f},2)-1
%         xvals2(1+(i-1)*xspacing:i*xspacing)=linspace(xSkelNew{f}(i),xSkelNew{f}(i+1),xspacing);
%     end
%     xvals2=xvals2(1:end-1);
%     if flipIDX(f)~=1
%         gradient=(ySkelNew{f}(flipIDX(f))-ySkelNew{f}(flipIDX(f)-1))/(xSkelNew{f}(flipIDX(f))-xSkelNew{f}(flipIDX(f)-1));
%     %     pp1=csape(xSkelNew{f}(1:flipIDX(f)),[0 ySkelNew{f}(1:flipIDX(f)) gradient],'clamped');
%     %     pp2=csape(xSkelNew{f}(flipIDX(f):end),[gradient ySkelNew{f}(flipIDX(f):end) 0],'clamped');
%         pp1=csape(xSkelNew{f}(1:flipIDX(f)),[ySkelNew{f}(1:flipIDX(f))],'variational');
%         pp2=csape(xSkelNew{f}(flipIDX(f):end),[ySkelNew{f}(flipIDX(f):end)],'variational');
%         v1=ppval(pp1,xvals1);
%         v2=ppval(pp2,xvals2);
%         hold on
% 
%         idX1=find(0<xvals1 & xvals1<width); xvals1=xvals1(idX1);
%         idX2=find(0<xvals2 & xvals2<width); xvals2=xvals2(idX2);
%         idv1=find(0<v1 & v1<height); v1=v1(idv1);
%         idv2=find(0<v2 & v2<height); v2=v2(idv2);
%         plot(xvals1,v1,'b',xvals2,v2,'g');
%         backbone{f}=[xvals1' v1'; xvals2' v2'];
%     else
%         pp=csape(xSkelNew{f}(flipIDX(f):end),[ySkelNew{f}(flipIDX(f):end)],'variational');
%         v=ppval(pp,xvals2);
%         hold on
% 
%         xvals2=xvals2(0<xvals2 & xvals2<width & 0<v & v<height);
%         v=v(0<v & v<height & 0<xvals2 & xvals2<width);
%         plot(xvals2,v,'b');
%         backbone{f}=[xvals2' v'];
%         hold off
%     end
%     
%     
% %     for i=1:Spacing:(numel(ySkelNew{f})-Spacing) 
% %         theta=atan2((ySkelNew{f}(i+Spacing)-ySkelNew{f}(i)),(xSkelNew{f}(i+Spacing)-xSkelNew{f}(i)))
% %         thetanew=atan2((ySkelNew{f}(i+1)-ySkelNew{f}(i)), (xSkelNew{f}(i+1)-xSkelNew{f}(i)));
% %         thetaend=atan2((ySkelNew{f}(i+Spacing)-ySkelNew{f}(i+Spacing-1)), (xSkelNew{f}(i+Spacing)-xSkelNew{f}(i+Spacing-1)));
% %         Rforward=[1 0 0; 0 cos(-theta) -sin(-theta); 0 sin(-theta) cos(-theta)];
% %         Rbackward=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
% %         for j=i:i+(Spacing)
% %             Xtemp{j}=Rforward*[0; xSkelNew{f}(j); ySkelNew{f}(j)];
% %             xtemp(j)=Xtemp{j}(2); ytemp(j)=Xtemp{j}(3);
% %         end
% %         if i==1
% %             startCond=0;
% %             endCond=thetaend;
% %         else
% %             startCond=tan(thetanew+oldtheta);
% %             endCond=(thetaend);
% %         end
% %             
% %         pp=csape(xtemp(i:i+(Spacing)), [ startCond ytemp(i:i+(Spacing)) endCond],'complete');
% %         vtemp=ppval(pp,linspace(xtemp(i),xtemp(i+(Spacing)),50));
% %         v=Rbackward*[zeros(length(vtemp),1) linspace(xtemp(i),xtemp(i+(Spacing)),50)' vtemp']';
% %         xnew=v(2,:); vnew=v(3,:);
% %         plot(xnew,vnew,'b');
% %         oldtheta=atan2((vnew(end)-vnew(end-1)),(xnew(end)-xnew(end-1)));
% %         oldCond=endCond;
% %     end
% %           
%     
%     Repeat=input(sprintf('Would you like to repeat the selection for frame %d? \n Yes (1) or no (2)',f),'s');
%     if Repeat=='1'
%         close
%     elseif Repeat=='2'
%         
%     else
%         error('Invalid input. Exiting...');
%     end
%     end
%     
%     % for each segment create normal lines
%     hold off
%     imshow((imageBW{f}));
%     hold on
%     ContX=ContnextX{1,f}./Pix-Offset;
%     ContX=width-ContX;
%     xl=xlim; yl=ylim;
%     axis manual
%     plot(ContX, ContnextY{1,f}./Pix-Offset,'r')
%     Contour{f}=[ContX ContnextY{1,f}./Pix-Offset];
% 
%     for i=1:(size(backbone{f},1)-1)
%         Grad{f}(i)=(backbone{f}(i+1,2)-backbone{f}(i,2))/(backbone{f}(i+1,1)-backbone{f}(i,1));
%         InvGrad{f}(i)=-1/(Grad{f}(i));
%         MidPoint{f}(i,:)=[(backbone{f}(i+1,1)+backbone{f}(i,1))/2 (backbone{f}(i+1,2)+backbone{f}(i,2))/2];
%         Intercept{f}(i)=MidPoint{f}(i,2)-Grad{f}(i)*MidPoint{f}(i,1);
% %         InvIntercept{f}(i)=(Grad{f}(i)-Intercept{f}(i))*MidPoint{f}(i,1)+Intercept{f}(i);
%         InvIntercept{f}(i)=MidPoint{f}(i,2)-InvGrad{f}(i)*MidPoint{f}(i,1);
% %         hold on
% %         refline(InvGrad{f}(i),InvIntercept{f}(i));
% %         hold off
%         Normals{f}(i,:)=[0 width (InvGrad{f}(i)*0+InvIntercept{f}(i)) (InvGrad{f}(i)*width+InvIntercept{f}(i))];
%         Intersections{f,i}=InterX([Normals{f}(i,1) Normals{f}(i,2); Normals{f}(i,3) Normals{f}(i,4)],[Contour{f}(:,1)'; Contour{f}(:,2)']);
%         IDX=knnsearch([Intersections{f,i}]',[MidPoint{f}(i,:)],'K',2);
%         Intersections{f,i}=Intersections{f,i}(:,IDX);
%         hold on
%         if isempty(Intersections{f,i})==0
%             line([Intersections{f,i}(1,1) Intersections{f,i}(1,2)],[Intersections{f,i}(2,1) Intersections{f,i}(2,2)]);
%             mesh{f}(i,:)=[Intersections{f,i}(1,1) Intersections{f,i}(2,1) Intersections{f,i}(1,2) Intersections{f,i}(2,2)];
%             DiametersSnake{f}(i,:)=sqrt((Intersections{f,i}(1,1)-Intersections{f,i}(1,2))^2 + (Intersections{f,i}(2,1)-Intersections{f,i}(2,2))^2);
%         end
%     end
%     DiametersSnaketd=std(DiametersSnake{f});
%     DiametersSnakeMean=mean(DiametersSnake{f});
%     idx=find((DiametersSnake{f}<(DiametersSnakeMean-2*DiametersSnaketd))|(DiametersSnake{f}>(DiametersSnakeMean+2*DiametersSnaketd)));
%     for i=1:numel(idx)
%         Intersections{f,idx(i)}=[];
%         DiametersSnake{f}(idx(i),:)=[];
%     end
%     
%     axis([xl yl]);
end
% close(h);
% save([name 'Mitomesh']);

set(0,'DefaultFigureWindowStyle','normal')