function testmesh(mesh,ContnextXsnake,ContnextYsnake)
    for i=1:size(mesh,2)
        figure
        hold on
        plot(ContnextXsnake{i},ContnextYsnake{i},'g')
        for j=1:size(mesh{1,i},1)
            line([mesh{1,i}(j,1) mesh{1,i}(j,3)],[mesh{1,i}(j,2) mesh{1,i}(j,4)])
        end
        axis equal
    end
end