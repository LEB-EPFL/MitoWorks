function choice=checkInitialMesh()
% Construct a questdlg with three options
    choice = questdlg('Do you already have the initial mesh?', ...
        'Initial mesh', ...
        'Yes - load initial mesh','No - run mitoMesh.m','Quit','Quit');
    % Handle response
    
end