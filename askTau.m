function answer =  askTau();

    prompt = {'Enter the smoothing factor for the backbone Tau:'};
    dlg_title = 'Enter Tau value';
    num_lines = 1;
    defaultans = {'10'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    
end