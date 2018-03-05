function IDXm=checkMesh(ContnextXsnake,ContnextYsnake,mesh,StartFrame,EndFrame)
IDXm=[];
set(0,'DefaultFigureWindowStyle','docked')

for f=StartFrame:EndFrame
           OK=0;
        while OK==0
            hold off
            plot(ContnextXsnake{f},ContnextYsnake{f},'g')
            hold on
            for j=1:size(mesh{1,f},1)
                line([mesh{1,f}(j,1) mesh{1,f}(j,3)],[mesh{1,f}(j,2) mesh{1,f}(j,4)])
            end           
            axis equal
            title(sprintf('Frame %d, correct (1) or wrong (2)?',f));
            choice=input(sprintf('Correct (1) or wrong(2) for frame %d?',f));
            
            if choice==1
                OK=1;
            elseif choice==2
                OK=1;
                IDXm=[IDXm f];
            else
            end
        end
end
set(0,'DefaultFigureWindowStyle','normal')
end
