function [Speed1, Speed2, Speed1proj, Speed2proj]= trackLeadingEdge(ContnextXsnakeDM1,ContnextYsnakeDM1, ContnextXsnakeDM2, ContnextYsnakeDM2,BreakFrame, EndFrame, MidPoint, N, Timestep, Direction, Pix);
% function measures the speed of retraction of the leading edge of daughter
% mitochondria after fission
%
%   INPUTS:
%   ContnextXsnakeDM1 and 2 = X coordinates of daughter mito contours
%   ContnextYsnakeDM1 and 2 = Y coordinates of daughter mito contours
%   BreakFrame = last frame before fission
%   EndFrame = last frame
%   MidPoint = last position of the constriction site before fission
%   N = number of nearest neighbour points to MidPoint defining the leading
%   edge for which the speed is measured
%   Timestep = time between frames in s
%   Direction = vector defining the orientation of the constriction site
%   Pix = pixel size in nm
% 
%   OUTPUTS:
%   Speed1 and 2 = speeds of leading edge retraction of daugher
%   mitochondria, with their corresponding X and Y components
%   Speed1proj and Speed2proj = projections of the speeds onto the
%   Direction vector


figure
axis equal

for f=BreakFrame:EndFrame
    hold off
    
    % find N nearest contour points
    [IDX1]=knnsearch([ContnextXsnakeDM1{f} ContnextYsnakeDM1{f}],MidPoint,'k',N);
    [IDX2]=knnsearch([ContnextXsnakeDM2{f} ContnextYsnakeDM2{f}],MidPoint,'k',N);
    Points1{f}=[ContnextXsnakeDM1{f}(IDX1) ContnextYsnakeDM1{f}(IDX1)];
    Points2{f}=[ContnextXsnakeDM2{f}(IDX2) ContnextYsnakeDM2{f}(IDX2)];
    
    plot(ContnextXsnakeDM1{f},ContnextYsnakeDM1{f},'g-', 'Linewidth',2)
    hold on
    plot(MidPoint(1),MidPoint(2),'kx','MarkerSize',18,'Linewidth',4)
    plot(ContnextXsnakeDM2{f},ContnextYsnakeDM2{f},'g-','Linewidth',2)
    plot(ContnextXsnakeDM1{f}(IDX1),ContnextYsnakeDM1{f}(IDX1),'b.','MarkerSize',18)
    plot(ContnextXsnakeDM2{f}(IDX2),ContnextYsnakeDM2{f}(IDX2),'b.','MarkerSize',18)
   
    
    % find their CM position
    CM1(f,:)=[mean(Points1{f}(:,1)) mean(Points1{f}(:,2))]'; % in pixels
    CM2(f,:)=[mean(Points2{f}(:,1)) mean(Points2{f}(:,2))]'; % in pixels
    
    plot(CM1(f,1),CM1(f,2),'rx','MarkerSize',18,'Linewidth',4)
    plot(CM2(f,1),CM2(f,2),'rx','MarkerSize',18,'Linewidth',4)
    
    Direction=normr(Direction); 
    
    drawnow;
    pause(0.25);
    
    if f~=BreakFrame
        Speed1(f,:)=[(CM1(f,1)-CM1(f-1,1)) (CM1(f,2)-CM1(f-1,2))]*Pix/Timestep;
        Speed2(f,:)=[(CM2(f,1)-CM2(f-1,1)) (CM2(f,2)-CM2(f-1,2))]*Pix/Timestep;
        
        % project along direction of constriction site
        Speed1proj(f,:)=dot(Speed1(f,:),Direction);
        Speed2proj(f,:)=dot(Speed2(f,:),Direction);
    end
end

end
