function MitoID=askMitoID()
    prompt = {'Enter the ID of the mitochondrion/mitochondria in INCREASING order:'};
    dlg_title = 'Enter mito ID(s)';
    num_lines = 1;
    defaultans = {'1'};
    MitoID = inputdlg(prompt,dlg_title,num_lines,defaultans);
end