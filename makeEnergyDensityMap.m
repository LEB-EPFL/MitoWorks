
function [normSlice]=makeEnergyDensityMap(IDmin, averageCurvatures, segmentRadii, sweep,  mesh,ContnextXsnake, ContnextYsnake, StartFrame, EndFrame);

%sweep is extent over which energy density to be calculated and displayed 

% find energy density per slice 
maxSlice=zeros(EndFrame,1);
for k=StartFrame:EndFrame
  EnvC{:,k}=(averageCurvatures{1,k}(:,1)+averageCurvatures{1,k}(:,2))./2;
  TubeC{:,k}=(1./(segmentRadii{1,k}));
  TotalC{:,k}=EnvC{:,k}+TubeC{:,k};
  ED{:,k}=(TotalC{:,k}).^2;
  
 %find max from whole movie (within segment of interest) and normalize
 %(normalize your color)
% for n=1:length(ED{:,k})

if (IDmin(k)-sweep)<=0
    roi{:,k}= ED{1, k}([1:IDmin(k)+sweep], 1);
    maxSlice(k,1)=max(roi{:,k});
elseif (IDmin(k)+sweep)>=size(ED{1, k},1)
    roi{:,k}= ED{1, k}([IDmin(k)-sweep:end], 1);
    maxSlice(k,1)=max(roi{:,k});
else
    roi{:,k}= ED{1, k}([IDmin(k)-sweep:IDmin(k)+sweep], 1);
    maxSlice(k,1)=max(roi{:,k});
end
%   maxMov=max(maxSlice);%max slice ED from full movie
%   normSlice{:,k}=ED{1,k}./maxMov;     
   
end

maxMov=max(maxSlice);
for f=StartFrame:EndFrame
    normSlice{:,f}=ED{1,f}./maxMov; 
end

  
%plotting

for i =StartFrame:EndFrame
   
   h=figure;
  
    
   set(gcf,'color','white');
    set(gca,'color','white');
    set(gca,'visible','off')
     set(0,'defaultfigurecolor','w')
%    plot(mesh{1,i}(:, 1:2:3), mesh{1,i}(:,2:2:4), 'Color', 'white', 'linewidth', 1.35)
plot(ContnextXsnake{1,i},ContnextYsnake{1,i}, 'Color', 'black', 'linewidth', 1.35)
   hold on
   if (IDmin(i)-sweep)<=0
       IDminrep=1;
   else
       IDminrep=(IDmin(i)-sweep);
   end
   
   if(IDmin(i)+sweep)>=size(mesh{1,i},1)
       IDmaxrep=size(mesh{1,i},1)-1;
   else
       IDmaxrep=(IDmin(i)+sweep-1);
   end
 for j=IDminrep:IDmaxrep
        set(gcf,'color','white');
    set(gca,'color','white');
       set(gca,'visible','off')
        set(0,'defaultfigurecolor','w')
     
    line([mesh{1,i}(j,1:2:3)],[mesh{1,i}(j,2:2:4)]);
    bounds_x=[mesh{1,i}(j,1); mesh{1,i}(j,3);mesh{1,i}(j+1, 3);mesh{1,i}(j+1,1) ];
    bounds_y=[mesh{1,i}(j,2); mesh{1,i}(j,4);mesh{1,i}(j+1, 4);mesh{1,i}(j+1,2) ];  
    
    C = normSlice{1,i}(j, 1);
    

 
    
    %fill each segment (plot color)
   set(0,'defaultfigurecolor','w')
   caxis manual
   caxis([0 1]);
   fill(bounds_x, bounds_y, C );
%    caxis manual
 
   colormap((jet));colorbar;
h=gcf;
h.InvertHardcopy = 'off';

% saveas(h, sprintf('frame%i.png', i));

    %fill entire roi

 end
   
end

end
