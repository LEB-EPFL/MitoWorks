function S=askSmoothingFactor();

    prompt = {'Enter the value of the smoothing factor for the contour generation (smoothing S*10 [nm]):'};
    dlg_title = 'Enter smoothing factor value';
    num_lines = 1;
    defaultans = {'17'};
    S = inputdlg(prompt,dlg_title,num_lines,defaultans);
end