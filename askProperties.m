function answer =  askProperties();

    prompt = {'Enter the start frame:', 'Enter the end frame', 'Enter the pixel size [nm]','Enter timestep between frames in [s]', 'Enter the frame before fission', 'Enter the temperature [deg C]'};
    dlg_title = 'Enter properties';
    num_lines = 1;
    defaultans = {'1','1','30','1', '1', '37'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    
end