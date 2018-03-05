function [k0 k1 k1bb]=generateCurvature(ContnextXsnake, ContnextYsnake, MidPoint, PixelSize,StartFrame,EndFrame)

%output k1 is the envelope curvature
h = waitbar(0,'Generating curvature values...');

for jj=StartFrame:EndFrame
        waitbar((jj-StartFrame)/(EndFrame-StartFrame));

  v{jj}=[ContnextXsnake{:,jj}, ContnextYsnake{:,jj}];
  Lines{jj}=[(1:(size(v{1,jj},1)-2))' (3:size(v{1,jj},1))'];% the curvature at a single point uses 3 neighbouring points (7 total: itself and 3 on each side)
    
  [k0{jj}]=LineCurvature2D(v{jj}, Lines{jj});
  [k1{jj}]=[k0{jj}./PixelSize];
  
  %backbone envelope curvature
  vB{jj}=[MidPoint{:,jj}(:,1),MidPoint{:,jj}(:,2)];
  LinesB{jj}=[(1:(size(vB{1,jj},1)-2))' (3:size(vB{1,jj},1))'];
  
  [k0bb{jj}]=LineCurvature2D(vB{jj}, LinesB{jj});
  [k1bb{jj}]=[k0bb{jj}./PixelSize];
    
end
close(h);
end

