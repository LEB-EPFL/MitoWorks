function Eta = askViscosity()

    prompt = {'Enter the viscosity value [kg m^{-1} s^{-1}]'};
    dlg_title = 'Enter viscosity value';
    num_lines = 1;
    defaultans = {'3e-2'};
    Eta = inputdlg(prompt,dlg_title,num_lines,defaultans);
end