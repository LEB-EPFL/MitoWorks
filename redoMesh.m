function [mesh, Contour, subNormals, DiametersSnake, Midpoint]=redoMesh(SegmentedPM,ContnextXsnake,ContnextYsnake,PixelSize, Offset,startFrame,endFrame, pos,meshD,IDXm,mesh, Contour, subNormals, DiametersSnake, Midpoint)

for i=1:length(IDXm)
    [meshtemp, Contourtemp, subNormalstemp, DiametersSnaketemp,MidPointtemp]=mitoMesh(SegmentedPM,ContnextXsnake,ContnextYsnake,PixelSize, Offset,IDXm(i),IDXm(i), pos,meshD);
    f=IDXm(i);
    mesh{1,f}=meshtemp{1,f};
    Contour{1,f}=Contourtemp{1,f};
    DiametersSnake{1,f}=DiametersSnaketemp{1,f};
    subNormals{1,f}=subNormalstemp{1,f};
    Midpoint{1,f}=MidPointtemp{1,f};
    close all;
end

end