function [FileName,PathName]=loadInitialMesh()
    [FileName,PathName] = uigetfile('*.mat','Load the initial mesh file');
end
