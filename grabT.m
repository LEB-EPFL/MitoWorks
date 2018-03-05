function Tres=grabT(Drp1,StartFrame, EndFrame, Timescale,FissionOrReversal);
% allows the user to select where Tf is for fission events and Tres is for
% reversals 
%
%   INPUTS:
%   Drp1 = Drp1 signal over time
%   StartFrame = first frame
%   EndFrame = last frame
%   Timescale = time between frames in seconds
%   FissionOrReversal = 1 if fission, 0 if reversal
%
%   OUTPUTS:
%   Tres = Tf for fission or [start finish] Tres for reversals

figure
plot((StartFrame-EndFrame:0)*Timescale,Drp1,'Linewidth',2)
xlabel('Time before fission [s]')
ylabel('Normalized Drp1 intensity')
title('Residency time')
if FissionOrReversal ==1
    [Tres,~]=ginput(1);
    Times=(StartFrame-EndFrame:0)*Timescale; Times=Times';
    ID=knnsearch(Times,Tres);
    Tres=Times(ID);
else
    [Tres,~]=ginput(2);
    Times=(StartFrame-EndFrame:0)*Timescale; Times=Times';
    ID=knnsearch(Times,Tres);
    Tres=Times(ID);
end

close;

end