function FramesToCorrect = askFrames();

x = inputdlg('Enter frames as space-separated numbers:',...
             'Frames to correct', [1 50]);
FramesToCorrect = str2num(x{:}); 
   

end

