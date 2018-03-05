function [tres, TresL]=grabT(Drp1,StartFrame, EndFrame, Timescale,FissionOrReversal, MinDiameter);
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

% 


figure
yyaxis left
% Drp1
plot(Drp1,'Linewidth',2); hold on; plot(Drp1, 'd');
xlabel('Time before fission [s]')
ylabel('Normalized Drp1 intensity')
title('Residency time')
yyaxis right
% MinDiameter
plot(MinDiameter,'Linewidth',2)
ylabel('MinDiameter')
% if FissionOrReversal ==1
%     [tres,~]=ginput(1);
% %     Times=(StartFrame-EndFrame:0)*Timescale; Times=Times';
% Times=StartFrame:1:EndFrame;Times=Times';
%     ID=knnsearch(Times,tres);
%     tres=Times(ID);
%     TresL=EndFrame-tres;
% else
%     [tres,~]=ginput(2);
%    % Times=(StartFrame-EndFrame:0)*Timescale; Times=Times';
%     Times=StartFrame:1:EndFrame;Times=Times';
% 
%     ID=knnsearch(Times,tres);
%     tres=Times(ID);
%     TresL=tres(2)-tres(1);
%     
% end
pause;
TresL=1;%TresL*Timescale;
tres=1'%tres*Timescale;

close;

end