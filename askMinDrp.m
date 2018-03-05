function Drp1FinalFrame = askMinDrp();

x = inputdlg('What is the final frame at which Drp1 is at the constriction site?',...
             'Drp1 final frame of interest', [1 30]);
Drp1FinalFrame = str2num(x{:}); 
   
end