function choice=checkFinalMesh()
% Construct a questdlg with three options
    choice = questdlg('Do you already have the final mesh?', ...
        'Initial mesh', ...
        'Yes - load final mesh','No - run mitoMesh.m','Quit','Quit');
    % Handle response
    
end