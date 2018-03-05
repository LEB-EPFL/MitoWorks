function [ TimeStamp, MinDiameter, IDmin] = MinDiameterSearch2(DiametersSnake, StartFrame, FissionFrame, TimeRes,Px,  ContnextXsnake, ContnextYsnake, subNormals, mesh)
%Output:TimeStamp is a matrix with times, where time = 0 is frame
%immediately before fission; MinDiameter is a matrix (same size as
%TimeStamp) where each row is the min diameter in nm, ind is temporary
%Input:
%IDX_fwhm,MinDiameter_fwhm (input)
%, DiametersFWHM,ContourFWHMside1 (output)

Total=FissionFrame-StartFrame+1; 
displays=ceil(Total/15); %15 subplots (frames) per display

   

for d=1:displays
 count=0;
    figure;
    
% MinD=zeros(FissionFrame,1);
% IDmin=zeros(FissionFrame,1);
% TimeStamp=zeros(FissionFrame,1);
% MinDiameter=zeros(FissionFrame,1);
for f=StartFrame+15*(d-1):StartFrame+15*d-1
    
    if f<=FissionFrame
          count=1+count;
%     [f count]
%       for count=1:15;
    [MinD(f,1), IDmin(f,1)]= min(DiametersSnake{1,f});
    
    TimeStamp(f,1)=f*TimeRes-((FissionFrame)*TimeRes); %frame immediately before fission is 0
    
%    figure;
%     axis equal
%     subplot(round(EndFrame-StartFrame),round((EndFrame-StartFrame)/2),(EndFrame-f+1))
%         subplot(round(EndFrame-StartFrame),round((EndFrame-StartFrame)),(EndFrame-f+1))
     

    subplot(3,5,count);

    plot(ContnextXsnake{1,f}(:,1), ContnextYsnake{1,f}(:,1), 'm-', 'linewidth', 2);
    hold on;
    line(subNormals{1,f}(IDmin(f,1),1:2:3),subNormals{1,f}(IDmin(f,1),2:2:4), 'LineWidth', 1.8, 'Color', [0 0 0])
    hold off;
    axis equal;
    drawnow;
    str = sprintf('frame %d ', f);
    title(str)
    
   end
%     end
% MinDiameter=Px.*MinD;
% IDmin;
    
    
    
end
    suptitle('Press any key once you are finished inspecting frames')

% MinDiameter=Px.*MinD;
% IDmin;

pause
% 
%ask user if min selection is OK
choice = questdlg('Is the minimum diameter location correct?', ...
    'Minimum Selection', ...
    'Yes, it is ok for all frames!','No, some frames need correction', 'Yes, it is ok for all frames!');
% Handle response
switch choice
    case 'Yes, it is ok for all frames!'
%         MinDiameter=zeros(length(StartFrame+15*(d-1):StartFrame+15*d-1),1);
        for f=StartFrame+15*(d-1):StartFrame+15*d-1
            if f<=FissionFrame
            MinDiameter(f,1)=Px.*MinD(f,1);
            IDmin(f,1); 
            end
        end
%         exp=1;
%             close all;
%             figure;
%             plot(TimeStamp(StartFrame:EndFrame,1),MinDiameter(StartFrame:EndFrame,1), 'k*')
%             title('Diameter versus time')
%             xlabel('time (s)')
%             ylabel('diameter (nm)')
%             
%         MinDiameter=MinDiameter(StartFrame:EndFrame,1);

       
    case 'No, some frames need correction'
         
       
%         exp=2; 
        FramesToCorrect=[];
          for f=StartFrame+15*(d-1):StartFrame+15*d-1
              if f<=FissionFrame
            MinDiameter(f,1)=Px.*MinD(f,1);
            IDmin(f,1); 
              end
          end
         FramesToCorrect=askFrames();
%     MinD_new=zeros(max(FramesToCorrect),1);
%     IDmin_new=zeros(max(FramesToCorrect),1);
         for j=1:length(FramesToCorrect)
             
             figure
           
             plot(ContnextXsnake{1,FramesToCorrect(j)},ContnextYsnake{1,FramesToCorrect(j)}, 'm-')
            title('Drag a square around region with min diameter')
             axis equal
             h=imrect;
             posSelect=getPosition(h);
%              close
             ind{:,j} = find(posSelect(1)<mesh{1,FramesToCorrect(j)}(:,1)&mesh{1,FramesToCorrect(j)}(:,1)<posSelect(1)+posSelect(3) & posSelect(2)<mesh{1,FramesToCorrect(j)}(:,2)&mesh{1,FramesToCorrect(j)}(:,2)<posSelect(2)+posSelect(4));

            [MinD_new(FramesToCorrect(j),1), IDmin_new(FramesToCorrect(j),1)]= min(DiametersSnake{1,FramesToCorrect(j)}(ind{:,j}));
%             
           close
            figure
            axis equal
            plot(ContnextXsnake{1,FramesToCorrect(j)}(:,1), ContnextYsnake{1,FramesToCorrect(j)}(:,1), 'r-', 'linewidth', 2);
            for i=1:length(ind{:,j})
            line(subNormals{1,FramesToCorrect(j)}(ind{:,j}(i),1:2:3),subNormals{1,FramesToCorrect(j)}(ind{:,j}(i),2:2:4), 'LineWidth', 0.7, 'Color', [0 0 0])
            end

            hold on
            
            line(subNormals{1,FramesToCorrect(j)}(ind{:,j}(IDmin_new(FramesToCorrect(j),1)),1:2:3),subNormals{1,FramesToCorrect(j)}(ind{:,j}(IDmin_new(FramesToCorrect(j),1)),2:2:4), 'LineWidth', 1.8, 'Color', [0 1 0])
             str = sprintf('Minimum selection corrected for frame %d, press any key to continue ', FramesToCorrect(j));
            title(str)
            axis equal
            drawnow;
            
            MinDiameter(FramesToCorrect(j),1)=Px.*MinD_new(FramesToCorrect(j),1);
            IDmin(FramesToCorrect(j), 1)=IDmin_new(FramesToCorrect(j),1)+min(ind{:,j})-1;
            close
         end
%          pause
%          figure;
%            plot(TimeStamp(StartFrame:EndFrame,1),MinDiameter(StartFrame:EndFrame,1), 'k*')
%            hold on
%             plot(TimeStamp(StartFrame:EndFrame,1),MinDiameter(StartFrame:EndFrame,1), 'k-')
%             title('Diameter versus time (snake contour)')
%             xlabel('time (s)')
%             ylabel('diameter (nm)')   
             
end
end




       
 


%find the min diameters from the FWHM contours
%knnsearch(X,Y) finds the nearest neighbor in X for each point in Y. 

% for t=StartFrame:EndFrame
% 
%     CandidatesIDX_fwhm{t}=knnsearch(ContourFWHMside1{1,t}(:,:),[mesh{1,t}(IDmin(t),1), mesh{1,t}(IDmin(t),2)],'k',10); %ContourFWHMside1 has same number of points as side 2 and width matrix
%     [minfwhm(t),idc(t)]=min(DiametersFWHM{1,t}(CandidatesIDX_fwhm{t}));
%     IDX_fwhm(t)= find(ContourFWHMside1{1,t}(:,:)==ContourFWHMside1{1,t}(CandidatesIDX_fwhm{t}(idc(t))));
%     MinDiameter_fwhm(t)= Px.*DiametersFWHM{1,t}(IDX_fwhm(t));
% end
%  figure;
%            plot(TimeStamp(StartFrame:EndFrame),MinDiameter_fwhm(StartFrame:EndFrame), 'b*')
%            hold on
%              plot(TimeStamp(StartFrame:EndFrame),MinDiameter_fwhm(StartFrame:EndFrame), 'b-')
%             title('Diameter versus time (FWHM)')
%             xlabel('time (s)')
%             ylabel('diameter (nm)')   
%  figure;
%            plot(TimeStamp(StartFrame:EndFrame),MinDiameter_fwhm(StartFrame:EndFrame), 'b*')
%            hold on
%            plot(TimeStamp(StartFrame:EndFrame,1),MinDiameter(StartFrame:EndFrame,1), 'k*')
%                   
%             title('Compare diameter versus time')
%             xlabel('time (s)')
%             ylabel('diameter (nm)') 
%             legend('FWHM', 'Snake Contour')
%             
close all;



plot(TimeStamp,MinDiameter);
xlabel('Time before fission [s]');
ylabel('Constriction diameter [nm]');


end