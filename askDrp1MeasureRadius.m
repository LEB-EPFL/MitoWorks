function Drp1Rad = askDrp1MeasureRadius();

x = inputdlg('What radius (nm) would you like to use to measure Drp1?',...
             'Drp1 measurement Radius', [1 50]);
Drp1Rad = str2num(x{:}); 
   

end

