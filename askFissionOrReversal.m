function FissionOrReversal=askFissionOrReversal()
    % Construct a questdlg with three options
    choice = questdlg('Is this a fission or a reversal event?', ...
        'Event type', ...
        'Fission','Reversal','Quit','Quit');
    % Handle response
    switch choice
        case 'Fission'
            FissionOrReversal = 1;
        case 'Reversal'
            FissionOrReversal = 0;
        case 'Quit'
            error('No mito found. Quitting...');
    end
end