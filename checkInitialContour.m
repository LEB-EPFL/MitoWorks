function choice=checkInitialContour()
% Construct a questdlg with three options
    choice = questdlg('Do you already have the initial contour?', ...
        'Initial contour', ...
        'Yes','No','Quit','Quit');
    % Handle response
    
end