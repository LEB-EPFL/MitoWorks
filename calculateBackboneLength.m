function length=calculateBackboneLength(MidPointDM2,BreakFramePlus1)
    length=0;
    for f=2:size(MidPointDM2{BreakFramePlus1},1)
        length=length+sqrt((MidPointDM2{BreakFramePlus1}(f,1)-MidPointDM2{BreakFramePlus1}((f-1),1))^2+(MidPointDM2{BreakFramePlus1}(f,2)-MidPointDM2{BreakFramePlus1}((f-1),2))^2);
    end
end
