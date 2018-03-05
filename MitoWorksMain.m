%function MitoWorkMain
tic;
addpath('TIFFStack'); addpath('code');
set(0,'defaultfigurecolor','white');

addpath(pwd);
set(0,'DefaultFigureWindowStyle','normal');

%% load segmented image
disp('Loading segmeted image...');
    [SegmentedFileName, SegmentedPathName]=loadSegmentedImage;
    [Seginfo, SegnumberOfFrames,Segwidth,Segheight] = getImageInfo([SegmentedPathName SegmentedFileName]);
%% load original image of mito
disp('Loading original image...');
    [OriginalSingleColorFileName, OriginalSingleColorPathName]=loadOriginalSingleColorImage;
    [Oriinfo, OrinumberOfFrames,Oriwidth,Oriheight] = getImageInfo([OriginalSingleColorPathName OriginalSingleColorFileName]);

    
%% ask user file save destination
    SaveFolder=uigetdir('', 'Save results in the selected folder');
    save([SaveFolder '\workspace']);

%% ask user for properties
disp('Inputing image properties...');
    Properties= askProperties();
    FissionOrReversal=askFissionOrReversal();
    save([SaveFolder '\workspace']);
%% track mito
disp('Identifying mito ID...');
    MitoFound=0;
    MitoFound=askMitoFound();

    
    while MitoFound==0
        disp('Tracking mitochondria using mitoTrack.m...');
        [Speeds, Tracks, segmentedMitoIDImageOld] = mitoTrack([SegmentedPathName SegmentedFileName],str2double(Properties{4}),false,1,str2double(Properties{1}), str2double(Properties{5}), str2double(Properties{3}));
        MitoFound=askMitoFound();
    end
    MitoID=askMitoID();
            close all;

    save([SaveFolder '\workspace']);

%% output specific mito ID mask
    disp('Creating segmented mask using mitoTrack.m...');
    [Speed, Track, segmentedMitoIDImageOld] = mitoTrack([SegmentedPathName SegmentedFileName],str2double(Properties{4}),true,str2double(MitoID{1}), str2double(Properties{1}), str2double(Properties{5}), str2double(Properties{3}));
    close all;
    save([SaveFolder '\workspace']);

%% apply mito mask to original segmented image
segmentedMitoIDImageOldDil=cell(1,str2double(Properties{5}));
OrigImage=cell(1,str2double(Properties{5}));
for f=str2double(Properties{1}):str2double(Properties{5})
    
    segmentedMitoIDImageOldDil{f}=imdilate(segmentedMitoIDImageOld{f},strel('ball',round(300/str2double(Properties{3})),0));
    OrigImage{f}=imread([OriginalSingleColorPathName OriginalSingleColorFileName],f);
end
    
   [SegmentedOriginalImage]=applyMaskToOriginalImage(segmentedMitoIDImageOldDil,[OriginalSingleColorPathName OriginalSingleColorFileName], str2double(Properties{1}),str2double(Properties{5}));
   [SegmentedPM]=applyMaskToOriginalImage(segmentedMitoIDImageOld,[SegmentedPathName SegmentedFileName], str2double(Properties{1}),str2double(Properties{5}));

    close all;
%% create contour
disp('Creating initial contour...');
% check if already exists
%     choice=checkInitialContour();
choice='No';
    switch choice
        case 'Yes'
            % if yes load the existing contour
            disp([choice ': initial contour exists. Loading initial contour...'])
            [InitialContourFileName, InitialContourPathName]=loadInitialContour();
            load([InitialContourPathName InitialContourFileName]);
        case 'No'
            % if no generate initial contour
            disp([choice ': initial contour does not exist. Generating initial contour using mitoCont.m...']);
            [ContnextXsnake, ContnextYsnake, pos]=generateSnakeContour(SegmentedPM,[OriginalSingleColorPathName OriginalSingleColorFileName],str2double(Properties{1}), str2double(Properties{5}),str2double(Properties{3}));
        case 'Quit'
            error('Quitting...');
    end
    %close all;
% if yes run it
% if no make it
    save([SaveFolder '\workspace']);

%% create first mesh
close all;
% check if already exists
disp('Creating initial mesh...');
% choice=checkInitialMesh();
choice='No - run mitoMesh.m';
switch choice
    % if yes load it
    case 'Yes - load initial mesh'
        disp([choice ': initial mesh exists. Loading initial mesh...'])
        [InitialMeshFileName,InitialMeshPathName]=loadInitialMesh();
        load([InitialMeshPathName InitialMeshFileName]);
        
    % if no make it
    case 'No - run mitoMesh.m'
        SpurCutoff=1;
        disp([choice ': initial mesh does not exist. Generating initial mesh using mitoMesh.m...']);
        [mesh, Contour, subNormals, DiametersSnake, MidPoint]=mitoMesh_auto(SegmentedPM,ContnextXsnake,ContnextYsnake,str2double(Properties{3}), 0,str2double(Properties{1}),str2double(Properties{5}), pos,3,SpurCutoff);
        IDXm=checkMesh(ContnextXsnake,ContnextYsnake,mesh,str2double(Properties{1}),str2double(Properties{5}));
        [mesh, Contour, subNormals, DiametersSnake, MidPoint]=redoMesh(SegmentedPM,ContnextXsnake,ContnextYsnake,str2double(Properties{3}), 0,str2double(Properties{1}),str2double(Properties{5}), pos,5,IDXm,mesh, Contour, subNormals, DiametersSnake, MidPoint);
    case 'Quit'
        error('Quitting...');
end

clear choice;
%     save([SaveFolder '\workspace']);

%% measure curvature along contour

origSubNormals=subNormals;
origMesh=mesh;
origDiametersSnake=DiametersSnake;

for f=str2double(Properties{1}):str2double(Properties{5})
        while (subNormals{1,f}(1,:)==[0 0 0 0])
%             display('cond');
            subNormals{1,f}(1,:)=[];
            mesh{1,f}(1,:)=[];
            DiametersSnake{1,f}(1)=[];
        end
end
    
%%
close all;
[k0, k1, k1bb]=generateCurvature(ContnextXsnake, ContnextYsnake,MidPoint, str2double(Properties{3}),str2double(Properties{1}),str2double(Properties{5}));    save([SaveFolder '\workspace']);
%     save([SaveFolder '\workspace']);

%% create FWHM profiles, measure FWHM diameters and calculate intensity corrected diameters 
% [ContourFWHMside1, ContourFWHMside2, DiametersFWHM, DiametersIntensity] = measureFWHMIntensityDiameter(SegmentedOriginalImage,pos,subNormals,mesh,1,ContnextXsnake, ContnextYsnake, str2double(Properties{1}), str2double(Properties{5}));
close all;  
save([SaveFolder '\workspace']);

% FWHM curvatures
%  [ContourFWHMside1, ContourFWHMside2, DiametersFWHM, DiametersIntensity, ConMito,dist,Int,fitG] = measurefitFWHMIntensityDiameter(SegmentedOriginalImage,pos,subNormals,mesh,1,ContnextXsnake, ContnextYsnake, str2double(Properties{1}), str2double(Properties{5}));


% [ks1, ks2]=generateCurvatureFWHM(ContourFWHMside1, ContourFWHMside2, str2double(Properties{3}),str2double(Properties{1}),str2double(Properties{5}));
%% Measure min diameter versus time
close all;
% [ TimeStamp, MinDiameter, IDmin] = MinDiameterSearch(DiametersSnake, StartFrame, FissionFrame, TimeRes,Px,  ContnextXsnake, ContnextYsnake, subNormals, mesh)
 [TimeStamp, MinDiameter, IDmin]=MinDiameterSearch2(DiametersSnake, str2double(Properties{1}), str2double(Properties{5}),str2double(Properties{4}),str2double(Properties{3}), ContnextXsnake, ContnextYsnake, subNormals, mesh);

        save([SaveFolder '\workspace']);

    
        %% calculate bending energy
close all;
Delta=10;
[energy, energyDensity, averageCurvatures, segmentRadii, IDmin, curvature, segCurv1, segCurv2, Extents, meshCurvature, energyInner, energyDensityInner, meshspacing, segLength, smallestID, subsegment, subsegmentArea, subsegmentInner, subsegmentAreaInner]=measureBendingEnergy(str2double(Properties{1}), str2double(Properties{5}), ContnextXsnake, ContnextYsnake, mesh, k1, IDmin, subNormals, 15,1, str2double(Properties{3}),SaveFolder,str2double(Properties{4}),Delta);
set(0,'defaultfigurecolor','white');
%     save([SaveFolder '\workspace']);

%% make a bending energy density movie
close all;
[normSlice]=makeEnergyDensityMap(IDmin, averageCurvatures, segmentRadii, 30,mesh,ContnextXsnake, ContnextYsnake, str2double(Properties{1}), str2double(Properties{5}));
%     save([SaveFolder '\workspace']);

%% extract minimum envelope curvatures ()
close all;
% [MinECconSite, smoothMinECconSite, avMinECconSite]=MinECurvatureExtract(segCurv1, segCurv2, str2double(Properties{1}),str2double(Properties{5}));
[TubeRnew, EnvRnew, EnvBackbone]=MinECurvatureExtract(segCurv1, segCurv2, str2double(Properties{1}),str2double(Properties{5}), segLength, smallestID, segmentRadii,Extents, k1bb);

%     save([SaveFolder '\workspace']);

%% FWHM contour + curvatures
close all;
set(0,'defaultfigurecolor','white');
[ContourFWHMside1, ContourFWHMside2, DiametersFWHM, DiametersIntensity, ConMito, dist, Int, fitG,ContourSegment,meshSegment, distO,c, MitoIntInt_av] = measurefitFWHMROIprofileSNR(SegmentedOriginalImage,pos,subNormals,mesh,1/3, ContnextXsnake, ContnextYsnake,Extents,50, str2double(Properties{1}), str2double(Properties{5}), segmentedMitoIDImageOld,MidPoint);
[ks1, ks2, ksbb,Side1, Side2, ContourFWHMccw, IDXside1, IDXside2]=generateCurvatureFWHM(ContourFWHMside1, ContourFWHMside2, str2double(Properties{3}),str2double(Properties{1}),str2double(Properties{5}));

    save([SaveFolder '\workspace']);

%%
close all; 
set(0,'defaultfigurecolor','white');

CandidatesIDX_fwhm=cell(1,str2double(Properties{5}));
minfwhm=zeros(1,str2double(Properties{5}));
idc=zeros(1,str2double(Properties{5}));
IDX_fwhm=zeros(1,str2double(Properties{5}));
MinDiameterFWHM=zeros(1,str2double(Properties{5}));
for t=str2double(Properties{1}):str2double(Properties{5})

    CandidatesIDX_fwhm{t}=knnsearch( ContourFWHMside1{1,t}(:,:),[mesh{1,t}(IDmin(t),1), mesh{1,t}(IDmin(t),2)],'k',10); %ContourFWHMside1 has same number of points as side 2 and width matrix
    [minfwhm(t),idc(t)]=min(DiametersFWHM{1,t}(CandidatesIDX_fwhm{t}));
    IDX_fwhm(t)= find(ContourFWHMside1{1,t}(:,:)==ContourFWHMside1{1,t}(CandidatesIDX_fwhm{t}(idc(t))));
    MinDiameterFWHM(t)= str2double(Properties{3}).*DiametersFWHM{1,t}(IDX_fwhm(t));
end
                                                                                                                                                                                                                                    % StartFrame, EndFrame, ContourFWHMside1, ContourFWHMside2, ks1,ks2, IDX_fwhm,PixelSize, SaveFolder,Delta
 [energyFWHM, energyDensityFWHM, averageCurvaturesFWHM, segmentRadiiFWHM, IDX_fwhm2, curvatureFWHM,ExtentFWHM, energyFWHMIn, energyDensityFWHMIn, segCurv1FWHM, segCurv2FWHM, meshspacingFWHM, segLengthFWHM, smallestIDFWHM, subsegmentFWHM, subsegmentAreaFWHM, subsegmentInnerFWHM, subsegmentAreaInnerFWHM]=measureBendingEnergyFWHM(str2double(Properties{1}), str2double(Properties{5}), ContourFWHMside1, ContourFWHMside2, ks1,ks2, IDX_fwhm, str2double(Properties{3}), SaveFolder, Delta);
                                                                                    % (segCurv1, segCurv2, StartFrame, EndFrame, segLength, smallestID, segmentRadii, Extents)
[TubeR_FWHM, EnvR_FWHM, EnvBackbone_FWHM]=MinECurvatureExtract(segCurv1FWHM, segCurv2FWHM, str2double(Properties{1}),str2double(Properties{5}), segLengthFWHM, smallestIDFWHM, segmentRadiiFWHM, ExtentFWHM, ksbb);

ConsLength=measureConsLength(averageCurvaturesFWHM,str2double(Properties{1}),str2double(Properties{5}),ContourFWHMside1, ContourFWHMside2, 7, IDX_fwhm, str2double(Properties{3}));

for i=str2double(Properties{1}):str2double(Properties{5})
        ContourSide1{i}=mesh{i}(:,1:2);
        ContourSide2{i}=mesh{i}(:,3:4);
    end

ConsLengthSnake=measureConsLength(averageCurvatures,str2double(Properties{1}),str2double(Properties{5}),ContourSide1, ContourSide2, 7, IDmin, str2double(Properties{3}));

% for i=1:size(ContourSegment,2)
% %     IDfwhm(i)=find(DiametersFWHM{i}==min(DiametersFWHM{i}));
%     IDfwhm(i)=round(size(ContourSegment{i},1)/2);
%     MinDiameterFWHM(i)=DiametersFWHM{i}(IDfwhm(i));
% end

%     save([SaveFolder '\workspace']);


%%  env curv FWHM
close all; 
% [MinECconSiteFWHM, smoothMinECconSiteFWHM, avMinECconSiteFWHM]=MinECurvatureExtract(segCurv1FWHM, segCurv2FWHM, str2double(Properties{1}),str2double(Properties{5}));

% R_Env1FWHM = 1/avMinECconSiteFWHM(end,1);
% R_Env2FWHM = 1/avMinECconSiteFWHM(end,2);
% R_EnvAvFWHM=(R_Env1FWHM+R_Env2FWHM)/2; %%these radii of curvature are in nm
% TubeEnvelopeFWHM=[R_Env1FWHM R_Env2FWHM R_EnvAvFWHM MinDiameterFWHM(end)];
%     save([SaveFolder '\workspace']);

%% bending energy for FWHM
% [energyFWHM, energyDensityFWHM, averageCurvaturesFWHM, segmentRadiiFWHM, IDX_fwhm, curvatureFWHM]=measureBendingEnergyFWHM(str2double(Properties{1}), str2double(Properties{5}), ContourFWHMside1, ContourFWHMside2, mesh, ks1,ks2, IDX_fwhm, subNormals, 15,1, str2double(Properties{3}), SaveFolder);
close all;
   save([SaveFolder '\workspace']);

% ======================================== DRP1 INTENSITY ==========================================================

%% Load Drp1 image
close all;
disp('Loading Drp1 image...');
[Drp1_FileName, Drp1_PathName]=loadDrpImage;
    save([SaveFolder '\workspace']);

%% Measure Drp1 flux at minimum constriction site

[CircleMask,Drp1IntensityConsSite,R0c, InMask] = measureDrpConsSite(ContnextXsnake, ContnextYsnake, IDmin, subNormals,TimeStamp, str2double(Properties{1}), str2double(Properties{5}), [Drp1_PathName Drp1_FileName], str2double(Properties{3}));
set(0,'defaultfigurecolor','white');

%subtract background and bleach correct Drp signal at constriction site
DrpIm=[Drp1_PathName Drp1_FileName];
MitoIm=[OriginalSingleColorPathName OriginalSingleColorFileName];
frameStamp_ConSite=str2double(Properties{1}):str2double(Properties{5});

%%
close all;
[Drp1IntensityConsSite_BleachBkgrndCorrected, BleachSlope,BleachYInter, BBox]= BackgroundBleachCorrect(str2double(Properties{1}), str2double(Properties{5}), DrpIm, MitoIm, Drp1IntensityConsSite, frameStamp_ConSite, CircleMask , R0c, InMask);
set(0,'defaultfigurecolor','white');
%     save([SaveFolder '\workspace']);

    Tres=grabT(Drp1IntensityConsSite_BleachBkgrndCorrected,str2double(Properties{1}), str2double(Properties{5}), str2double(Properties{4}),FissionOrReversal);

%% Load Drp1 Tracks and select tracks within frame range
disp('Loading Drp1 tracks...');
[Drp1Tracks_FileName, Drp1Tracks_PathName]=importTracks;
[Drp1Tracks, RadiusTrack] =  FilterTracks([Drp1Tracks_PathName Drp1Tracks_FileName], str2double(Properties{1}), str2double(Properties{5}));
%     save([SaveFolder '\workspace']);

%% user selects tracks of interest by drawing a polygon 
[SelectedTracks] = SelectTracks([OriginalSingleColorPathName OriginalSingleColorFileName],[Drp1_PathName Drp1_FileName],str2double(Properties{1}), str2double(Properties{5}), Drp1Tracks, str2double(Properties{3}) );
[SelectedTracksSpeeds] = trackSelectedTracks(SelectedTracks,str2double(Properties{4}));
%     save([SaveFolder '\workspace']);

%% Measure Drp1 flux from final Drp1 track (final Drp1 spot that leads to fission)

% find the final Drp1 track of interest
disp('Final Drp1 frame...');
choice=questdlg('Is Drp1 visible at the final frame you selected?', ...
        'Final Drp1 frame', ...
        'Yes','No -insert final frame','No -insert final frame');
switch choice  
    case 'Yes'
        [DrpTrackCons]=IsolateDrpCons(SelectedTracks,  str2double(Properties{5}), IDmin, subNormals, str2double(Properties{3}));
       [SpeedDrpTrackCons]=trackFinalTrack(DrpTrackCons,str2double(Properties{4}));


    case 'No -insert final frame'
       
Drp1FinalFrame = askMinDrp();

[DrpTrackCons]=IsolateDrpCons(SelectedTracks,  Drp1FinalFrame, IDmin, subNormals, str2double(Properties{3}));
[SpeedDrpTrackCons]=trackFinalTrack(DrpTrackCons,str2double(Properties{4}));
end

% measure the Drp1 integrated intensity of the final track
disp('Drp1 measurement radius');
Drp1Rad=askDrp1MeasureRadius;
[Drp1IntensityFinalTrack] = measureDrpSignal_finalTrack([Drp1_PathName Drp1_FileName], DrpTrackCons, Drp1Rad, str2double(Properties{3}), TimeStamp,  str2double(Properties{4}), str2double(Properties{1}), str2double(Properties{5}) );

%bleach correct the final track
frameStamp_final=DrpTrackCons(1,4):DrpTrackCons(end,4);
% DrpImage=[Drp1_PathName Drp1_FileName];
% MitoImage=[OriginalSingleColorPathName OriginalSingleColorFileName];
[Drp1IFinalTrack_BleachCorrected]= LinearBleachCorrection(Drp1IntensityFinalTrack, frameStamp_final,BleachSlope,BleachYInter);


% measure the Drp1 integrated intensity of multiple Drp1 tracks that enter
% the ROI, within their mask
[Drp1Mask, Drp1IntensityMultiTrack] = measureDrpSignal_multiTrack([Drp1_PathName Drp1_FileName],SelectedTracks,Drp1Rad,str2double(Properties{3}),str2double(Properties{1}),str2double(Properties{5}), CircleMask);

%bleach correct the merged track
[frameStamp_multiTracks] = frameStamper(SelectedTracks,str2double(Properties{5}));
[Drp1ImultiTracks_BleachCorrected]= LinearBleachCorrection(Drp1IntensityMultiTrack, frameStamp_multiTracks, BleachSlope,BleachYInter);
close all; 
save([SaveFolder '\workspace']);


%% curvature at Drp1

[EnvCurvDrp1_FinalTrack, EnvCurv_SelectedTracks,  MinEnvCurvatures ] = measureEnvCurvAtDrp1(str2double(Properties{1}), str2double(Properties{5}), ContnextXsnake, ContnextYsnake,DrpTrackCons, SelectedTracks, curvature, str2double(Properties{3}));
close all;    
% save([SaveFolder '\workspace']);
    
[EnvCurvDrp1_FinalTrackfwhm, EnvCurv_SelectedTracksfwhm,  MinEnvCurvaturesfwhm ] = measureEnvCurvAtDrp1FWHM(str2double(Properties{1}), str2double(Properties{5}),ContourFWHMside1, ContourFWHMside2,DrpTrackCons, SelectedTracks, curvatureFWHM, str2double(Properties{3}));
close all; 
%   save([SaveFolder '\workspace']);
% measure the envelope curvature around the final Drp1 track and the
% previously selected Drp1 tracks
% [EnvCurvDrp1_FinalTrack, EnvCurv_SelectedTracks,  MinEnvCurvatures] = measureEnvC(StartFrame, EndFrame, ContnextXsnake, ContnextYsnake,DrpTrackCons, SelectedTracks, curvature);

%     save([SaveFolder '\workspace']);

%% TO DO:Track Drp1 inside MitoWorksMain
% addpath('C:\Program Files\Fiji.app\scripts');
% Miji(false);
close all; 
% tfmin=nan;
% tfmax=nan;
% [param]=sigm_fit((-length(MinDiameter):-1)*str2double(Properties{4}),(Drp1IntensityConsSite_BleachBkgrndCorrected-min(Drp1IntensityConsSite_BleachBkgrndCorrected))/(max(Drp1IntensityConsSite_BleachBkgrndCorrected)-min(Drp1IntensityConsSite_BleachBkgrndCorrected)),[tfmin; tfmax; nan; nan]);
% tf=param(3);

x=(-length(MinDiameter(str2double(Properties{1}):str2double(Properties{5}))):-1)*str2double(Properties{4});
y=(Drp1IntensityConsSite_BleachBkgrndCorrected(str2double(Properties{1}):str2double(Properties{5}))-min(Drp1IntensityConsSite_BleachBkgrndCorrected(str2double(Properties{1}):str2double(Properties{5}))))/(max(Drp1IntensityConsSite_BleachBkgrndCorrected(str2double(Properties{1}):str2double(Properties{5})))-min(Drp1IntensityConsSite_BleachBkgrndCorrected(str2double(Properties{1}):str2double(Properties{5}))));
psig=fitSigmoid(x',y'); t12sig=psig(3);
pexp=fitShiftedExp(x',y'); t12exp=(log(0.5)/pexp(2))-1;

% save([SaveFolder '\workspace']);

[LineScanI1, LineScanI2] = lineScan([OriginalSingleColorPathName OriginalSingleColorFileName],[Drp1_PathName Drp1_FileName],MidPoint, str2double(Properties{1}),str2double(Properties{5}),str2double(Properties{3}));
if FissionOrReversal==1
DrpPostFis = trackDrpPostFis(averageCurvaturesFWHM,[OriginalSingleColorPathName OriginalSingleColorFileName], [Drp1_PathName Drp1_FileName], str2double(Properties{5}), [ContourFWHMside1{str2double(Properties{5})} ContourFWHMside2{str2double(Properties{5})}],DrpTrackCons, str2double(Properties{3}))
end

%% 
close all;

figure
[axX, ~, ~]=plotyy((-length(MinDiameter):-1)*str2double(Properties{4}),MinDiameter,(-length(MinDiameter):-1)*str2double(Properties{4}),(Drp1IntensityConsSite_BleachBkgrndCorrected-min(Drp1IntensityConsSite_BleachBkgrndCorrected))/(max(Drp1IntensityConsSite_BleachBkgrndCorrected)-min(Drp1IntensityConsSite_BleachBkgrndCorrected)));
xlabel('Time before fission [s]');
ylabel(axX(1),'Constriction diameter [nm] (snake)');
ylabel(axX(2),'Normalized Drp1 intensity');
title('Method 1: cons');
savefig([SaveFolder '\DvsDrp1cons']);

figure
[axX, ~, ~]=plotyy((-length(MinDiameter):-1)*str2double(Properties{4}),MinDiameter,(-length(MinDiameter):-1)*str2double(Properties{4}),(Drp1ImultiTracks_BleachCorrected-min(Drp1ImultiTracks_BleachCorrected))/(max(Drp1ImultiTracks_BleachCorrected)-min(Drp1ImultiTracks_BleachCorrected)));
xlabel('Time before fission [s]');
ylabel(axX(1),'Constriction diameter [nm] (snake)');
ylabel(axX(2),'Normalized Drp1 intensity');
title('Method 2: multitrack');
savefig([SaveFolder '\DvsDrp1multitrack']);

% figure
% [axX axY1 axY2]=plotyy((-length(MinDiameter_fwhm):-1)*str2double(Properties{4}),MinDiameter_fwhm,(-length(MinDiameter_fwhm):-1)*str2double(Properties{4}),(Drp1IntensityConsSite_BleachBkgrndCorrected-min(Drp1IntensityConsSite_BleachBkgrndCorrected))/(max(Drp1IntensityConsSite_BleachBkgrndCorrected)-min(Drp1IntensityConsSite_BleachBkgrndCorrected)));
% xlabel('Time before fission [s]');
% ylabel(axX(1),'Constriction diameter [nm] (FWHM)');
% ylabel(axX(2),'Normalized Drp1 intensity');
% title('Method 1: cons');
% savefig([SaveFolder '\DFWHMvsDrp1cons']);
% 
% figure
% [axX axY1 axY2]=plotyy((-length(MinDiameter_fwhm):-1)*str2double(Properties{4}),MinDiameter_fwhm,(-length(MinDiameter_fwhm):-1)*str2double(Properties{4}),(Drp1ImultiTracks_BleachCorrected-min(Drp1ImultiTracks_BleachCorrected))/(max(Drp1ImultiTracks_BleachCorrected)-min(Drp1ImultiTracks_BleachCorrected)));
% xlabel('Time before fission [s]');
% ylabel(axX(1),'Constriction diameter [nm] (FWHM)');
% ylabel(axX(2),'Normalized Drp1 intensity');
% title('Method 2: multitrack');
% savefig([SaveFolder '\DFWHMvsDrp1multitrack']);
% 
% close all;

%% 
close all;
%    save([SaveFolder '\workspace']);
%% ======================================================== PULLING FORCE ==============================================

if (str2double(Properties{2}) > str2double(Properties{5})) && FissionOrReversal==1
close all;
%% track daughter mitochondria from fission frame
set(0,'defaultfigurecolor','white');

disp('Identifying daughter mito IDs...');
    MitoFound=0;
    MitoFound=askMitoFound();

    
    while MitoFound==0
        close all;
        disp('Tracking daughter mitochondria using mitoTrack.m...');
        [Speeds2, Tracks2,SegmentedMitoIDImage] = mitoTrack([SegmentedPathName SegmentedFileName],str2double(Properties{4}),false,1,str2double(Properties{5})+1, str2double(Properties{2}), str2double(Properties{3}));
        close all; clear Speeds2; clear Tracks2;
        MitoFound=askMitoFound();
    end
    MitoIDs=askMitoID();
    MitoIDs=str2num(MitoIDs{1,1});
%         save([SaveFolder '\workspace']);

%% track daughter mitochondria

% need to varify quantities

    disp('Creating merged segmented mask using mitoTrack.m...');
    [SpeedDM1, TrackDM1, SegmentedMitoIDImageDM1] = mitoTrack([SegmentedPathName SegmentedFileName],str2double(Properties{4}),true,(MitoIDs(1)), str2double(Properties{5})+1, str2double(Properties{2}),str2double(Properties{3}));
    [SpeedDM2, TrackDM2, SegmentedMitoIDImageDM2] = mitoTrack([SegmentedPathName SegmentedFileName],str2double(Properties{4}),true,(MitoIDs(2)), str2double(Properties{5})+1, str2double(Properties{2}), str2double(Properties{3}));
    close all;
    
    BreakingFrame=splitLastFrame(segmentedMitoIDImageOld,str2double(Properties{5}),IDmin,subNormals);
    
    PullingFrames = assignPullingFrames(BreakingFrame,SegmentedMitoIDImageDM1,SegmentedMitoIDImageDM2,str2double(Properties{5}),str2double(Properties{2}));
    [PullingFrames1, PullingFrames2]=splitMito(BreakingFrame,str2double(Properties{5}),str2double(Properties{2}),SegmentedMitoIDImageDM1,SegmentedMitoIDImageDM2);
%     save([SaveFolder '\workspace']);
%         save([SaveFolder '\workspace']);

    %% create individual images
    [SegmentedOriginalImageDM1]=applyMaskToOriginalImage(PullingFrames1,[SegmentedPathName SegmentedFileName], str2double(Properties{5}),str2double(Properties{2}));
    [SegmentedOriginalImageDM2]=applyMaskToOriginalImage(PullingFrames2,[SegmentedPathName SegmentedFileName], str2double(Properties{5}),str2double(Properties{2}));
%         save([SaveFolder '\workspace']);

    %% generate daughter mito contours
    [ContnextXsnakeDM1, ContnextYsnakeDM1, posDM1]=generateSnakeContour(SegmentedOriginalImageDM1,[OriginalSingleColorPathName OriginalSingleColorFileName],str2double(Properties{5}), str2double(Properties{2}),str2double(Properties{3}));
    close all;
    [ContnextXsnakeDM2, ContnextYsnakeDM2, posDM2]=generateSnakeContour(SegmentedOriginalImageDM2,[OriginalSingleColorPathName OriginalSingleColorFileName],str2double(Properties{5}), str2double(Properties{2}),str2double(Properties{3}));
    close all;
%         save([SaveFolder '\workspace']);

    %% create daughter mito mesh
    SpurCutoff=5;
        [meshDM1, ContourDM1, subNormalsDM1, DiametersSnakeDM1, MidPointDM1]=mitoMesh_auto(PullingFrames1,ContnextXsnakeDM1,ContnextYsnakeDM1,str2double(Properties{3}), 0,str2double(Properties{5}),str2double(Properties{2}), pos,3,SpurCutoff);
        IDXm=checkMesh(ContnextXsnakeDM1,ContnextYsnakeDM1,meshDM1,str2double(Properties{5}),str2double(Properties{2}));
        [meshDM1, ContourDM1, subNormalsDM1, DiametersSnakeDM1, MidPointDM1]=redoMesh(PullingFrames1,ContnextXsnakeDM1,ContnextYsnakeDM1,str2double(Properties{3}), 0,str2double(Properties{5}),str2double(Properties{2}), pos,3,IDXm,meshDM1, ContourDM1, subNormalsDM1, DiametersSnakeDM1, MidPointDM1);
    close all;
        [meshDM2, ContourDM2, subNormalsDM2, DiametersSnakeDM2, MidPointDM2]=mitoMesh_auto(PullingFrames2,ContnextXsnakeDM2,ContnextYsnakeDM2,str2double(Properties{3}), 0,str2double(Properties{5}),str2double(Properties{2}), pos,3,SpurCutoff);
        IDXm=checkMesh(ContnextXsnakeDM2,ContnextYsnakeDM2,meshDM2,str2double(Properties{5}),str2double(Properties{2}));
        [meshDM2, ContourDM2, subNormalsDM2, DiametersSnakeDM2, MidPointDM2]=redoMesh(PullingFrames2,ContnextXsnakeDM2,ContnextYsnakeDM2,str2double(Properties{3}), 0,str2double(Properties{5}),str2double(Properties{2}), pos,3,IDXm,meshDM2, ContourDM2, subNormalsDM2, DiametersSnakeDM2, MidPointDM2);    close all;
%         save([SaveFolder '\workspace']);    

%% track leading edge
    P1=[(mesh{str2double(Properties{5})}(IDmin(str2double(Properties{5}))-1,1)+mesh{str2double(Properties{5})}(IDmin(str2double(Properties{5}))-1,3))/2 (mesh{str2double(Properties{5})}(IDmin(str2double(Properties{5}))-1,2)+mesh{str2double(Properties{5})}(IDmin(str2double(Properties{5}))-1,4))/2];
    P2=[(mesh{str2double(Properties{5})}(IDmin(str2double(Properties{5}))+1,1)+mesh{str2double(Properties{5})}(IDmin(str2double(Properties{5}))+1,3))/2 (mesh{str2double(Properties{5})}(IDmin(str2double(Properties{5}))+1,2)+mesh{str2double(Properties{5})}(IDmin(str2double(Properties{5}))+1,4))/2];
    N=50; % number of nearest neighbours 
    [Speed1, Speed2, Speed1proj, Speed2proj]= trackLeadingEdge(ContnextXsnakeDM1,ContnextYsnakeDM1, ContnextXsnakeDM2, ContnextYsnakeDM2,str2double(Properties{5}), str2double(Properties{2}), [(mesh{str2double(Properties{5})}(IDmin(str2double(Properties{5})),1)+mesh{str2double(Properties{5})}(IDmin(str2double(Properties{5})),1))/2 (mesh{str2double(Properties{5})}(IDmin(str2double(Properties{5})),2) + mesh{str2double(Properties{5})}(IDmin(str2double(Properties{5})),4))/2],N,str2double(Properties{4}),[P2(:,1)-P1(:,1) P2(:,2)-P1(:,2)], str2double(Properties{3}));

    
    %% calculate pulling force
    close all;
    [DaughterSpeed, DaughterSpeedV]=mitoPull(PullingFrames,str2double(Properties{4}),str2double(Properties{1}),str2double(Properties{2}),str2double(Properties{5}), str2double(Properties{3}),[P2(:,1)-P1(:,1) P2(:,2)-P1(:,2)]);
    close all;
    
    
    L1=calculateBackboneLength(MidPointDM1,str2double(Properties{5})+1); % in pixels
        L2=calculateBackboneLength(MidPointDM2,str2double(Properties{5})+1); % in pixels

        D1=mean(DiametersSnakeDM1{str2double(Properties{5})+1}); % in pixels
        D2=mean(DiametersSnakeDM2{str2double(Properties{5})+1}); % in pixels
    

        Eta = askViscosity();
        Eta=str2double(Eta{1,1});
        [F1_VEM, F2_VEM, T_VEM] = viscoElasticModel(Speed1proj, Speed2proj, Eta, MidPointDM1, MidPointDM2, DiametersSnakeDM1, DiametersSnakeDM2, str2double(Properties{5}),str2double(Properties{2}), str2double(Properties{3}), MinDiameter);
        
        
        
        
        imshow(SegmentedOriginalImageDM1{1,str2double(Properties{5})+1});
        FricCoeff1 = calculateFricCoeff(Eta,D1*str2double(Properties{3})*1e-9,L1*str2double(Properties{3})*1e-9);
        close all;
        imshow(SegmentedOriginalImageDM2{1,str2double(Properties{5})+1});
        FricCoeff2 = calculateFricCoeff(Eta,D2*str2double(Properties{3})*1e-9,L2*str2double(Properties{3})*1e-9);
        close all;
    FricCoeff=[FricCoeff1 FricCoeff2];
    [DaughterForce, ForceError] = calculateForce(FricCoeff,DaughterSpeed,str2double(Properties{6})); % should be in [N]
    [DaughterForceX, ForceErrorX]=calculateForce(FricCoeff,[DaughterSpeedV{1}(:,1) DaughterSpeedV{2}(:,1)], str2double(Properties{6}));
    [DaughterForceY, ForceErrorY]=calculateForce(FricCoeff,[DaughterSpeedV{1}(:,2) DaughterSpeedV{2}(:,2)], str2double(Properties{6}));
    
    Tension=max(max(DaughterForce((str2double(Properties{5})+1):(str2double(Properties{5})+3),:)))/(min(nonzeros(MinDiameter))*1e-9);
    close all;    save([SaveFolder '\workspace']);

    %% making sure it's the same as in binFissions
    if (str2num(Properties{5})+1+round(3/str2double(Properties{4})))<str2num(Properties{2})
            maxPullingForceFnew=max(max(DaughterForce((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),:)));
            sumPullingForceFnew=sum(max(DaughterForce((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),:)));

            stdPullingForceFnew=std([DaughterForce((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),1); DaughterForce((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),2)]);
            maxPullingForceSUMFnew=max((DaughterForce((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),1)+DaughterForce((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),2)));
            stdPullingForceSUMFnew=std(DaughterForce((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),1)+DaughterForce((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),2));
            MaxSpeedLEproj=max(max(abs([Speed1proj((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4})))) Speed2proj((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))))])));
            SumSpeedLEproj=sum(max(abs([Speed1proj((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4})))) Speed2proj((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))))])));
            
            MaxSpeedLE=max(max(abs([sqrt(Speed1((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),1).^2+Speed1((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),2).^2) sqrt(Speed2((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),1).^2 + Speed2((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),2).^2)])));
            SumSpeedLE=sum(max(abs([sqrt(Speed1((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),1).^2+Speed1((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),2).^2) sqrt(Speed2((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),1).^2 + Speed2((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),2).^2)])));
%             CompSpeed(FNnew)=max([DaughterForce]
            CompSpeed=[abs(Speed1proj-DaughterSpeed(:,1)) abs(Speed2proj-DaughterSpeed(:,2))];
            MaxCompSpeedLE=max(max(CompSpeed((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),:)));
            stdMaxCompSpeedLE=std([CompSpeed((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),1); CompSpeed((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),2)]);
            sumCompSpeedLE=sum(max(CompSpeed((str2num(Properties{5})+1):(str2num(Properties{5})+1+round(3/str2double(Properties{4}))),:)));

            

        else
            maxPullingForceFnew=max(max(DaughterForce((str2num(Properties{5})+1):end,:)));
            sumPullingForceFnew=sum(max(DaughterForce((str2num(Properties{5})+1):end,:)));
            stdPullingForceFnew=std([DaughterForce((str2num(Properties{5})+1):end,1); DaughterForce((str2num(Properties{5})+1):end,2)]);
            maxPullingForceSUMFnew=max((DaughterForce((str2num(Properties{5})+1):end,1)+DaughterForce((str2num(Properties{5})+1):end,2)));
            stdPullingForceSUMFnew=std(DaughterForce((str2num(Properties{5})+1):end,1)+DaughterForce((str2num(Properties{5})+1):end,2));
            MaxSpeedLEproj=max(max(abs([Speed1proj((str2num(Properties{5})+1):end) Speed2proj((str2num(Properties{5})+1):(end))])));
            SumSpeedLEproj=sum(max(abs([Speed1proj((str2num(Properties{5})+1):end) Speed2proj((str2num(Properties{5})+1):(end))])));

            MaxSpeedLE=max(max(abs([sqrt(Speed1((str2num(Properties{5})+1):end,1).^2+Speed1((str2num(Properties{5})+1):end,2).^2) sqrt(Speed2((str2num(Properties{5})+1):(end),1).^2+Speed2((str2num(Properties{5})+1):(end),2).^2)])));
            SumSpeedLE=sum(max(abs([sqrt(Speed1((str2num(Properties{5})+1):end,1).^2+Speed1((str2num(Properties{5})+1):end,2).^2) sqrt(Speed2((str2num(Properties{5})+1):(end),1).^2+Speed2((str2num(Properties{5})+1):(end),2).^2)])));
            CompSpeed=[(Speed1proj-DaughterSpeed(:,1)) (Speed2proj-DaughterSpeed(:,2))];
            MaxCompSpeedLE=max(max(CompSpeed((str2num(Properties{5})+1):end,:)));
            stdMaxCompSpeedLE=std([CompSpeed((str2num(Properties{5})+1):end,1); CompSpeed((str2num(Properties{5})+1):end,2)]);
            sumCompSpeedLE=sum(max(CompSpeed((str2num(Properties{5})+1):end,:)));
    end
        
    
    TensionMW=Tension;
        TensionFnew=maxPullingForceFnew/(pi()*MinDiameter(end)*1e-9);
        SumTensionFnew=sumPullingForceFnew/(pi()*MinDiameter(end)*1e-9);

        stdTensionFnew=stdPullingForceFnew/(pi()*MinDiameter(end)*1e-9);
        ForceFperAnew=maxPullingForceFnew/subsegmentArea(str2double(Properties{5}),IDmin(str2double(Properties{5})));
        if (str2double(Properties{5})+1+round(3/str2double(Properties{4})))<=length(T_VEM)
            TensionF_VEMnew=max(T_VEM((str2double(Properties{5})+1):(str2double(Properties{5})+1+round(3/str2double(Properties{4})))));
            ForceF_VEMnew=max(max([F1_VEM((str2double(Properties{5})+1):(str2double(Properties{5})+1+round(3/str2double(Properties{4})))) F2_VEM((str2double(Properties{5})+1):(str2double(Properties{5})+1+round(3/str2double(Properties{4}))))]));
        else
            TensionF_VEMnew=max(T_VEM((str2double(Properties{5})+1):end));
            ForceF_VEMnew=max(max([F1_VEM((str2double(Properties{5})+1):end) F2_VEM((str2double(Properties{5})+1):end)]));
        end

        ForceF_VEMperA=ForceF_VEMnew/subsegmentArea(str2double(Properties{5}),IDmin(str2double(Properties{5})));

    
end
save([SaveFolder '\workspace']);

toc;
