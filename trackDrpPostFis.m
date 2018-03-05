function DrpPostFis = trackDrpPostFis(averageCurvatures,MitoImage, DrpImage, EndFrame, meshseg,Drp1track, Pixel)
% determines position of Drp1 post fission
%
%   INPUTS
%   averageCurvatures = [env tube] curvature cell for each frame
%   MitoImage = path to green channel
%   DrpImage = path to red channel
%   EndFrame = frame just before fission
%   meshseg = mesh used to create segments
%   Drp1track = positions of Drp1 punctae for each frame [x y]
%   Pixel = pixel size in nm
%
%   OUTPUTS:
%   DrpPostFis = string with Drp1 status
    
    % find max ED segment before fission
    f=EndFrame;
    ED=(averageCurvatures{f}(:,1)+averageCurvatures{f}(:,2)).^2;
    ID_EDmin=find(ED==max(ED));
    P(f,:)=[mean([meshseg(ID_EDmin+1,1) meshseg(ID_EDmin+1,3) meshseg(ID_EDmin,1) meshseg(ID_EDmin,3)]) mean([meshseg(ID_EDmin+1,2) meshseg(ID_EDmin+1,4) meshseg(ID_EDmin,2) meshseg(ID_EDmin,4)])];
    IDdrp=find(Drp1track(:,4)==EndFrame);
    
    % display image with previous positions of 
    figure; imshow(imfuse(imread(MitoImage,EndFrame+1),imread(DrpImage,EndFrame+1)),[])
    hold on
    title('o = Drp, x = max ED');
    plot(P(EndFrame,1),P(EndFrame,2),'kx','Markersize',18,'Linewidth',2)
    if isempty(IDdrp)==0
        plot(Drp1track(IDdrp,2)/Pixel,Drp1track(IDdrp,3)/Pixel,'ro','Linewidth',2);
    else
        warning('Warning: no Drp1 puncta detected.')
    end
    pause
    str={'Good side'; 'Bad side'; 'Both sides'; 'Neither side'; 'Don''t know'};
    [s,v] = listdlg('PromptString','Select an option:',...
                    'SelectionMode','single',...
                    'ListString',str);
    DrpPostFis=str{s};
    close;
end