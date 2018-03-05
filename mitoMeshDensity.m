function [mesh Contour subNormals DiametersSnake MidPoint] = mitoMeshDensity(imesh, StartFrame, EndFrame, Resize, ExtensionDistance,height,width, Contour);

clear iMidPoint iMeshSlope MidPointX MidPointY MeshSlope MidPoints SmoothBackbone

for f=StartFrame:EndFrame
    for j=1:size(imesh{f},1)
        iMidPoint{f,j}=[(imesh{f}(j,1)+imesh{f}(j,3))/2 (imesh{f}(j,2)+imesh{f}(j,4))/2];
%         iMeshSlope{f,j}=(imesh{f}(j,4)-imesh{f}(j,2))/(imesh{f}(j,3)-imesh{f}(j,1));
    end
    
    
    for j=1:(size(imesh{f},1)-1)
        for i=1:Resize
            if i==1
                MidPointX{f}(j,1)=(iMidPoint{f,j}(1));
                MidPointY{f}(j,1)=(iMidPoint{f,j}(2));
%                 MeshSlope{f}(j,1)=iMeshSlope{f,j};
            else
                MidPointX{f}(j,i)=(iMidPoint{f,j}(1)+iMidPoint{f,j+1}(1))*(i-1)/Resize;
                MidPointY{f}(j,i)=(iMidPoint{f,j}(2)+iMidPoint{f,j+1}(2))*(i-1)/Resize;
%                 MeshSlope{f}(j,i)=(iMeshSlope{f,j}+iMeshSlope{f,j+1})*(i-1)/Resize;
            end
        end
    end
    
    MidPointX{f}=reshape(MidPointX{f}',prod(size(MidPointX{f})),[]);
    MidPointY{f}=reshape(MidPointY{f}',prod(size(MidPointY{f})),[]);
%     MeshSlope{f}=reshape(MeshSlope{f}',prod(size(MeshSlope{f})),[]);
    MidPoints{f}=[MidPointX{f} MidPointY{f}];
    
    SmoothBackbone{f}=lssmooth([MidPoints{f}(:,1) MidPoints{f}(:,2)],50);
    figure
    hold on
    plot(iMidPoint{f}(:,1), iMidPoint{f}(:,2),'g-', 'LineWidth',2);
%     plot(MidPoints{f}(:,1), MidPoints{f}(:,2), 'r-', 'LineWidth',2);
%     plot(SmoothBackbone{f}(:,1),SmoothBackbone{f}(:,2),'b-', 'LineWidth',2);
    NoSegments=size(SmoothBackbone{f},1)-1;

    
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

set(0,'DefaultFigureWindowStyle','normal');

