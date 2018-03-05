function I=findCnstArea(subsegmentArea,Time,Area);
    
% initialize
    i=1;
    AreaP=subsegmentArea(Time,i);
    AreaC=subsegmentArea(Time,i+1)+AreaP;
    DiffC=abs(AreaC-Area);
    DiffP=abs(AreaP-Area);
    if AreaC<Area
    % repeat until error starts increasing
        while DiffC<DiffP
            i=i+1
            AreaP=AreaP+subsegmentArea(Time,i);
            AreaC=AreaC+subsegmentArea(Time,i+1);
            DiffC=abs(AreaC-Area);
            DiffP=abs(AreaP-Area);
        end

       I=i; 
    else
        I=1;
    end
    
end
