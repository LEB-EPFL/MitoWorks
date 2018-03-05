function [LineScanI1, LineScanI2] = lineScan(MitoImage,DrpImage,MidPoint, StartFrame,EndFrame,Pixel)
% calculates the intensities of the mito and drp channel along the backbone
% of the mesh
%
%   INPUTS
%   MitoImage = path to green channel
%   DrpImage = path to red channel
%   MidPoint = each cell {f} for frame f contains the [x y] coordinates of
%               backbone (produced from mesh)
%   StartFrame = first frame 
%   EndFrame = end frame
%   Pixel = pixel size in nm
%
%   OUTPUTS
%   LinescanI1 and 2 = [I x y] measured intensity from image along linescan
%                       at points with coordinates x and y in nm

for f=StartFrame:EndFrame
    % read image
    Imito=imread(MitoImage,f);
    Idrp=imread(DrpImage,f);
    
    % measure intensity along backbone
    [cmito_x cmito_y cmito_new] = improfile(Imito,MidPoint{f}(:,1),MidPoint{f}(:,2));
    [cdrp_x cdrp_y cdrp_new] = improfile(Idrp,MidPoint{f}(:,1),MidPoint{f}(:,2));
    LineScanI1{f}(:,1)=cmito_new;
    LineScanI2{f}(:,1)=cdrp_new;
    LineScanI1{f}(:,2)=cmito_x*Pixel; LineScanI1{f}(:,3)=cmito_y*Pixel;
    LineScanI2{f}(:,2)=cdrp_x*Pixel; LineScanI2{f}(:,3)=cdrp_y*Pixel;
end

end