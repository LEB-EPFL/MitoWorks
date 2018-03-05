function Rad_curvDetection= askRadiusCurvatureDetection()

x = inputdlg('Specify cut-off radius (nm) for envelope curvature measurement',...
             'Envelope curvature measurement Radius', [1 50]);
         
Rad_curvDetection = str2num(x{:}); 
   
end

