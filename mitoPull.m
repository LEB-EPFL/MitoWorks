function [Speed, SpeedV]=mitoPull(InputImage,timestep,StartFrame,EndFrame,BreakFrame, pixel_size, ConstrictionLine)
% splits the parent mitochondrion in two at the constriction site and
% calculates the recoil speeds of daughter mitochondria post fission

% Inputs: 
% InputImage = segmented image of parent + daughter mito around the fission
%              frame
% timestep = in [s]
% StartFrame
% EndFrame
% BreakFrame = frame before fission occurs
% pixel_size = in [nm]
% ConstrictionLine = line running along the constriction site
%
% Outputs: 
% Speed = in [nm] of daughter mitochondria
imageBW=cell(1,EndFrame);
CC=cell(1,EndFrame);
imageBW_labeled=cell(1,EndFrame);
MitoMeasurements=cell(1,EndFrame);
NumberOfMito=cell(1,EndFrame);
for k = BreakFrame:EndFrame
    imageBW{k} = InputImage{k};
    [CC{k}, imageBW_labeled{k}] = labelBWImage(imageBW{k},0);
    MitoMeasurements{k} = regionprops(imageBW_labeled{k}, imageBW{k}, 'all');    
    NumberOfMito{k}=size(MitoMeasurements{k},1);
    if k==BreakFrame
        MaxNumberOfMito=NumberOfMito{k};
    else
        if NumberOfMito{k}>MaxNumberOfMito;
            MaxNumberOfMito=NumberOfMito{k};
        end
    end
end
fission=0;

h = waitbar(0,'Running mitoPull.m');
Speed=zeros(EndFrame,MaxNumberOfMito);
SpeedV=cell(1,MaxNumberOfMito);
for Frame=BreakFrame:EndFrame
        waitbar((Frame-BreakFrame)/(EndFrame-BreakFrame));

    if (NumberOfMito{Frame}==1)
%                     imshow(imageBW{k});
%                     drawnow;
%                     pause(0.5);

        if (fission==1)
            break;
            error('A mito has left the FOV.');
        end
    elseif (NumberOfMito{Frame}==2)
        fission=1;
        if (Frame==BreakFrame)  
            %perpendicular to diameter line
            ConstrictionLine=normr(ConstrictionLine);
          hold on
            labels{1}='1'; labels{2}='2';
%             figure
            imshow(insertObjectAnnotation(imageBW_labeled{BreakFrame},'rectangle',[MitoMeasurements{Frame}(1).BoundingBox; MitoMeasurements{Frame}(2).BoundingBox],labels,'TextBoxOpacity',0.9,'FontSize',18));
            drawnow;
%             pause(1);
            hold off
        else
            % find closest mito
            [mitoID]=knnsearch([MitoMeasurements{Frame-1}(1).Centroid; MitoMeasurements{Frame-1}(2).Centroid],[MitoMeasurements{Frame}(1).Centroid; MitoMeasurements{Frame}(2).Centroid]);
            tempMeasurements=MitoMeasurements{Frame};
            MitoMeasurements{Frame}(1)=tempMeasurements(mitoID(1));
            MitoMeasurements{Frame}(2)=tempMeasurements(mitoID(2));
            labels{1}='1'; labels{2}='2';
%             figure
            hold on
            imshow(insertObjectAnnotation(imageBW_labeled{Frame},'rectangle',[MitoMeasurements{Frame}(1).BoundingBox; MitoMeasurements{Frame}(2).BoundingBox],labels,'TextBoxOpacity',0.9,'FontSize',18));
            drawnow;
%             pause(1);
        hold off
            for Mito=1:NumberOfMito{Frame} 
                Speed(Frame,Mito)=dot(((MitoMeasurements{Frame}(Mito).Centroid-MitoMeasurements{Frame-1}(Mito).Centroid)/timestep),ConstrictionLine)*pixel_size;
                SpeedV{Mito}(Frame,:)=(MitoMeasurements{Frame}(Mito).Centroid-MitoMeasurements{Frame-1}(Mito).Centroid)*pixel_size/timestep;
            end
        end
    else
%         Frame
        error('More than 2 mito found');
    end
end
close(h);
% cd(name)
% close
% figure
% plot(1:length(indX),(Speed(indX,1)),'ro-',1:length(indX),Speed(indX,2),'bo-')
% xlabel('Frame post fission');
% ylabel('Speed [nm/s]');
% grid on
% savefig([name '_SpeedV2']);
% 
% figure
% plot([1:length(nonzeros(Force(:,1)))],(nonzeros(Force(:,1))-min(nonzeros(Force(:,1))))/(max(nonzeros(Force(:,1)))-min(nonzeros(Force(:,1)))),'ro-',[1:length(nonzeros(Force(:,2)))],nonzeros(Force(:,2))/N,'bo-')
% xlabel('Frame post fission');
% ylabel('Force [a. u.]');
% grid on
% savefig([name '_ForceV2']);


% fprintf(output_file,'% %12s %12s %12s &12s %12s %12s %12s %12s \r\n','Speed1','Speed2','Force1','Force2','SpeedMitoBefore','Timestep','PixSize','Bcoeff');
% Y=[[1:noFrames]' Speed(:,1) Speed(:,2) Force(:,1) Force(:,2) SpeedMitoBefore ones(noFrames,1)*timestep ones(noFrames,1)*pixel_size ones(noFrames,1)*Bcoeff];
% fprintf(output_file,'%f %.4e %.4e %.4e %.4e %.4e %f %f %f \r\n',Y');
% fclose(output_file);
% save([name '_data']);
% cd ..;