function [EnvCurvDrp1_FinalTrack, EnvCurv_SelectedTracks,  MinEnvCurvatures ] = measureEnvCurvAtDrp1(StartFrame, EndFrame, ContnextXsnake, ContnextYsnake,DrpTrackCons, SelectedTracks, curvature, Pixel)
%measureEnvC determines the envelope curvature around the Drp1 spots using
%Drp1 tracks

Rad_curvDetection=askRadiusCurvatureDetection(); 
Rad=Rad_curvDetection./Pixel;

%for all selected tracks
EnvCurv_SelectedTracks=cell(length(SelectedTracks),EndFrame);
frames=cell(length(SelectedTracks),1);
LastFrame=zeros(length(SelectedTracks),1);
% newLastFrame
xyPOS=cell(length(SelectedTracks),1);

for l=1:length(SelectedTracks)%for each track
      frames{l,1}=SelectedTracks{l,1}(:,4); %each track has a matrix containing number of all frames
          LastFrame(l,1)=SelectedTracks{l,1}(size(SelectedTracks{l,1},1),4); %last frame for each track

      if LastFrame(l,1)>EndFrame
          newLastFrame{l,:}=knnsearch(SelectedTracks{l,1}(:,4), EndFrame, 'k', 2);
          if isempty(newLastFrame{l,:}==EndFrame)==1
          LastFrame(l,1)=min(newLastFrame(l,1));
          else
           LastFrame(l,1)=EndFrame;
          end
      end
      
      xyPOS{l,1}=SelectedTracks{l,1}(:,2:3)./Pixel;  %positions of the Drp1 spots for each track in a single matrix; detection radius drawn around this

      for f=frames{l,1}(1,1):LastFrame(l,1)%for each frame
         delta=LastFrame(l,1)-frames{l,1}(1,1);
          %define the contour 
          Contour{1}= ContnextXsnake{1,f};
          Contour{2}= ContnextYsnake{1,f};
          for count=1:delta+1
          xy{f}= xyPOS{l,1}(count, :);
          end
         ContPointDrp1{l,f}=pointsincircle(Contour, Rad,xy{f});
  
          %plot circles and points of interest on contour
          figure;
          str=sprintf('Track %d, Frame %d',l,f);
        
            plot(ContnextXsnake{1,f},ContnextYsnake{1,f}, 'k-');
            hold on 
            plot(xy{f}(:,1),xy{f}(:,2), 'c+');
            hold on
            circlePLOT(xy{f}(:,1),xy{f}(:,2),Rad);
            hold on
            plot(ContPointDrp1{l,f}.in{1,1},ContPointDrp1{l,f}.in{1,2}, 'r.')
            title(str);
      end
end
       
%for each track, find the corresponding envelope curvatures on side 1
%and/or 2 of the contour
% IN_side1=cell(length(SelectedTracks),LastFrame(s,1));
% IN_side2=cell(length(SelectedTracks),LastFrame(s,1));
% idx1=cell(length(SelectedTracks),LastFrame(s,1));
% idx2=cell(length(SelectedTracks),LastFrame(s,1));
    for s=1:length(SelectedTracks)
    for fr=frames{s,1}(1,1):LastFrame(s,1)
   
%ismember(A,B) returns an array containing 1 (true) where the data in A is
%found in B, A is the contour(named 'curvature' here) and B is the matrix of detected points (ContPointDrp1) 

 A=curvature{1,fr}(:,1:2); %side 1 contour
 a=curvature{2,fr}(:,1:2); %side 2 contour

B=[ContPointDrp1{s,fr}.in{1,1} ContPointDrp1{s,fr}.in{1,2}];

[IN_side1{s,fr}]=ismember(A,B, 'rows'); 
[IN_side2{s,fr}]=ismember(a,B, 'rows');


[idx1{s,fr}(:,1), idx1{s,fr}(:,2)]=find(IN_side1{s,fr}==1);

[idx2{s,fr}(:,1), idx2{s,fr}(:,2)]=find(IN_side2{s,fr}==1);

EnvCurv_SelectedTracks{s,fr}{1}=curvature{1,fr}(idx1{s,fr}(:,1),:); %side 1, a cell array where s is the track, fr is the frame and contains the x,y and curvature locations of side {1}
EnvCurv_SelectedTracks{s,fr}{2}=curvature{2,fr}(idx2{s,fr}(:,1),:); %side 2 (same as above for side 2)
    end
    end
    
%for Drp1 at constriction site

DrpTrackID=DrpTrackCons(1,1);
target1=cell(length(SelectedTracks),1);
for m=1:length(SelectedTracks)
target1{m,1}=find(SelectedTracks{m,1}(1,1)==DrpTrackID);
end
target2=cellfun(@isempty, target1); %inverts target1 and moves to matrix format

targetRow=find(target2==0);
EnvCurvDrp1_FinalTrack=cell(1,EndFrame);
for n=StartFrame:EndFrame
    EnvCurvDrp1_FinalTrack{1,n}=EnvCurv_SelectedTracks{targetRow, n};
% n
end


NoTracks=length(SelectedTracks);
Mside1=cell(NoTracks,EndFrame);
RowI_1=cell(NoTracks,EndFrame);
yn1=cell(NoTracks,EndFrame);
ROI_1=cell(NoTracks,EndFrame);
IDXs1=cell(NoTracks,EndFrame);
MinEnvCurvatures=cell(NoTracks,EndFrame);
Mside2=cell(NoTracks,EndFrame);
RowI_2=cell(NoTracks,EndFrame);
yn2=cell(NoTracks,EndFrame);
ROI_2=cell(NoTracks,EndFrame);
IDXs2=cell(NoTracks,EndFrame);
for t=1:NoTracks
    for i= StartFrame:EndFrame
%        
        %side1
        if isempty(EnvCurv_SelectedTracks{t,i})==0  
            if isempty(EnvCurv_SelectedTracks{t,i}{1,1})==0
            [Mside1{t,i},RowI_1{t,i}]=min(EnvCurv_SelectedTracks{t,i}{1,1}(:,3)); %min curvature on side1
            
            %ROI_1 (row of interest) is min curvature entry (RowI) in curvature Side 1 matrix
            [yn1{t,i},ROI_1{t,i}]=ismember(EnvCurv_SelectedTracks{t,i}{1,1}(RowI_1{t,i}, :), curvature{1,i}, 'rows'); 
            
            %find nearest neaighbours to min curvature point
            IDXs1{t,i}=knnsearch(curvature{1,i}(:, 1:2), curvature{1,i}(ROI_1{t,i}, 1:2), 'k', 5);
            
            MinEnvCurvatures{t,i}{1,1}=(1/5)*sum(curvature{1,i}(IDXs1{t,i}, 3));
            
            end
        end
        
        %side2
        if isempty(EnvCurv_SelectedTracks{t,i})==0
            if length(EnvCurv_SelectedTracks{t,i}{1,2})~=0
                [Mside2{t,i},RowI_2{t,i}]=min(EnvCurv_SelectedTracks{t,i}{1,2}(:,3));%min curvature on side2
                
            %ROI_2(row of interest) is min curvature entry (RowI) in curvature Side 1 matrix
            [yn2{t,i},ROI_2{t,i}]=ismember(EnvCurv_SelectedTracks{t,i}{1,2}(RowI_2{t,i}, :), curvature{2,i}, 'rows'); 
            
            %find nearest neaighbours to min curvature point
            IDXs2{t,i}=knnsearch(curvature{2,i}(:, 1:2), curvature{2,i}(ROI_2{t,i}, 1:2), 'k', 5);
            
            MinEnvCurvatures{t,i}{1,2}=(1/5)*sum(curvature{2,i}(IDXs2{t,i}, 3));
            end
        end
         
        end
end
   
%plot min curvature for each side of each track
TrackSide1=cell(NoTracks,EndFrame);
TrackSide2=cell(NoTracks,EndFrame);
for g=1:NoTracks
    for b=StartFrame:EndFrame
        
        %side1
        if isempty(MinEnvCurvatures{g,b})==0
            if isempty(MinEnvCurvatures{g,b}{1,1})==0
                TrackSide1{g,b}=MinEnvCurvatures{g,b}{1,1};
%                empty1=cellfun('isempty', TrackSide1);
%                TrackSide1(empty1)={NaN}
            end
        end
        
        %side2
         if isempty(MinEnvCurvatures{g,b})==0
            if length(MinEnvCurvatures{g,b})==2
                TrackSide2{g,b}=MinEnvCurvatures{g,b}{1,2};
%                empty2=cellfun('isempty', TrackSide2);
%                TrackSide2(empty2)={NaN}
            end
        end
        
    
    end
end
% % TO DO: plot env curvatures 

% % for j=1:NoTracks
% %    
% %     for c=StartFrame:EndFrame
% %         if isempty(TrackSide1{j,c})==0
% % 
% %        TS1(j,c)=TrackSide1{j,c};
% %         countFrames(j,c)=c;
% %         
% %         end
% %    
% %     end
% %      figure
% %         plot(countFrames, TS1, 'x')
% %     
% % end
    


        
        
        


      

            

end



   
    
    
    




