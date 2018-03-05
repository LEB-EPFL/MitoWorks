function MitoFound=askMitoFound()
    % Construct a questdlg with three options
    choice = questdlg('Have you found the mitochondrion you were looking for?', ...
        'Finding mito', ...
        'Yes','No','Quit','Quit');
    % Handle response
    switch choice
        case 'Yes'
            disp([choice ': mito identified. Please identify the mito ID.'])
            MitoFound = 1;
        case 'No'
            disp([choice ': mito not found. Running MitoTrack.m ...'])
            MitoFound = 0;
        case 'Quit'
            error('No mito found. Quitting...');
    end
end