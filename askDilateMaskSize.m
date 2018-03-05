function answer = askDilateMaskSize();
    prompt = {'Enter by how much you would like to dilate the mask [pixels] (try 7 as default):'};
    dlg_title = 'Enter dilating mask size';
    num_lines = 1;
    defaultans = {'7'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
end
