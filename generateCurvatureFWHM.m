function [ks1 ks2, ksbb, Side1, Side2, ContourFWHMccw, IDXside1, IDXside2, ContourFWHM]=generateCurvatureFWHM(ContourFWHMside1, ContourFWHMside2, PixelSize,StartFrame,EndFrame)

%output k1 is the envelope curvature
h = waitbar(0,'Generating curvature values...');
%side1

% for f=StartFrame:EndFrame
%     ContourFWHM{1,f}=[ContourFWHMside1{1,f};flipud(ContourFWHMside2{1,f})];
% end

for jj=StartFrame:EndFrame
        waitbar((jj-StartFrame)/(EndFrame-StartFrame));
        
   ContourFWHM{1,jj}=[ContourFWHMside1{1,jj};flipud(ContourFWHMside2{1,jj})];
   [ContourFWHMccw{1,jj}(:,1),ContourFWHMccw{1,jj}(:,2)] =poly2ccw(ContourFWHM{1,jj}(:,1),ContourFWHM{1,jj}(:,2));
   
   
%    FrameLen(jj)=length(ContourFWHMside1{1,jj});

   v{jj}=[ContourFWHMccw{1,jj}(:,1), ContourFWHMccw{1,jj}(:,2)];
   Lines{jj}=[(1:(size(v{1,jj},1)-2))' (3:size(v{1,jj},1))'];% the curvature at a single point uses 3 neighbouring points (7 total: itself and 3 on each side)
    
  [k01{jj}]=LineCurvature2D(v{jj}, Lines{jj});
  [k02{jj}]=[k01{jj}./PixelSize]; %these are the env curvatures of the merged, ccw FWHM contour
  
  IDXside1{jj}=knnsearch(ContourFWHMccw{1,jj}, ContourFWHMside1{1,jj});
  IDXside2{jj}=knnsearch(ContourFWHMccw{1,jj}, ContourFWHMside2{1,jj});
  
  ks1{jj}=k02{jj}(IDXside1{jj});
  ks2{jj}=k02{jj}(IDXside2{jj});
  
%   ks1{jj}=k02{jj}(1:FrameLen(jj));
%   ks2{jj}=k02{jj}(FrameLen(jj)+1:end);
    
end
% close(h);

%side2
% for ii=StartFrame:EndFrame
%         waitbar((ii-StartFrame)/(EndFrame-StartFrame));
% 
%    v{ii}=[ContourFWHMside2{1,ii}(:,1), ContourFWHMside2{1,ii}(:,2)];
%    Lines{ii}=[(1:(size(v{1,ii},1)-2))' (3:size(v{1,ii},1))'];% the curvature at a single point uses 3 neighbouring points (7 total: itself and 3 on each side)
%     
%   [k02{ii}]=LineCurvature2D(v{ii}, Lines{ii});
%   [ks2{ii}]=[k02{ii}./PixelSize];
%   
% end

%asign envelope curvatures of k02 to each side




 %backbone envelope curvature for FWHM contour
 for t=StartFrame:EndFrame
  vB{t}=[0.5*(ContourFWHMside1{1,t}(:,1)+ContourFWHMside2{1,t}(:,1)), 0.5*(ContourFWHMside1{1,t}(:,2)+ContourFWHMside2{1,t}(:,2))];
  LinesB{t}=[(1:(size(vB{1,t},1)-2))' (3:size(vB{1,t},1))'];
  
  [ksbb1{t}]=LineCurvature2D(vB{t}, LinesB{t});
  [ksbb{t}]=[ksbb1{t}./PixelSize];
 end
close(h);

for f=StartFrame:EndFrame
    
    Side1{1,f}(:,1)=ContourFWHMccw{1,f}(IDXside1{1,f},1); 
    Side1{1,f}(:,2)=ContourFWHMccw{1,f}(IDXside1{1,f},2);
    
    Side2{1,f}(:,1)=ContourFWHMccw{1,f}(IDXside2{1,f},1); 
    Side2{1,f}(:,2)=ContourFWHMccw{1,f}(IDXside2{1,f},2);
end

end

