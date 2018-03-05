function [width, xLeft, xRight,yLeft,yRight,posLeftMax,posRightMax]= getFWHM(varargin)
% Function [width, xLeft, xRight,yLeft,yRight,posLeftMax,posRightMax]=
% getFWXM(x,y, 'options')
% finds full width at X* the max (X=0.5 for full width at  half max), 
% for (x,y) data than can have multiple peaks.
% options:
% - 'Half': Instead of FWHM, find FW at X*max
% - 'nbreaks': number of breakpoints for splinefitting to find
% maxima.
% - 'searchStep': step size for search for FWXM
% - 'FigZ' 


% Input: x,y,searchstep, half
narginchk(2,10);
x=varargin{1,1};
y=varargin{1,2};
nbreaks=10;
Half=0.5;
searchStep=0.1;
% loop to check and asign input variables
for i=3:2:size(varargin,2)-1
    switch varargin{1,i}
        case 'Half'
            Half=varargin{1,i+1};
        case 'nbreaks'
            nbreaks=varargin{1,i+1};
        case 'searchStep'
            searchStep=varargin{1,i+1};
        case 'FigZ'
            fig=varargin{1,i+1};
        otherwise
            error(['invalid input' varargin{1,i} ])
    end
end


% pp=splinefit(x,y,nbreaks); % splinefit has the bad property of creating artifacts at the border, esp. when y=0 for multiple consecutive points.
% spline=ppval(pp,x);

[maxI,maxX]=findpeaks(y,x,'MinPeakHeight',max(y)/2);

% from found peaks, filter out the real left and right peaks
switch length(maxX)
    case 0 % what if no peaks are found?
        % check if there is a significant edge max.
        
        % else
        posLeftMax=NaN;
        valLeftMax=NaN;
        posRightMax=NaN;
        valRightMax=NaN;
    case 1 % only one peak
        posLeftMax=maxX;
        valLeftMax=maxI;
        posRightMax=maxX;
        valRightMax=maxI;
    case 2 % two peaks
        if maxX(1)<maxX(2) % if maxX(1) is left peak
            posLeftMax=maxX(1);
            valLeftMax=maxI(1);
            posRightMax=maxX(2);
            valRightMax=maxI(2);
        else
            posLeftMax=maxX(2);
            valLeftMax=maxI(2);
            posRightMax=maxX(1);
            valRightMax=maxI(1);
        end
    otherwise % more than two peaks...
%         disp('getFWXM says: "More than 2 peaks, hard to pick two." ')
%         %too spammy
        % sort peaks, pick most left and most right peaks.
        [maxX, idx]=sort(maxX,2,'ascend');
        maxI=maxI(idx);
        posLeftMax=maxX(1);
        valLeftMax=maxI(1);
        posRightMax=maxX(end);
        valRightMax=maxI(end);
%         % Debug
%         figure,
%         plot(x,y)
%         hold on, plot(x,spline)
%         plot(posLeftMax,valLeftMax,'+')
%         plot(posRightMax,valRightMax,'+')
%         hold off
end

% From left peak, scan left to find half max
posLeft = posLeftMax;
valCur = Inf;
while posLeft>=x(1) && valCur > valLeftMax*Half
    posLeft =posLeft-searchStep;
    valCur = interp1(x,y,posLeft);
end
% From right peak, scan right to find half max
posRight = posRightMax;
valCur = Inf;
while posRight<=x(end) && valCur > valRightMax*Half
    posRight =posRight+searchStep;
    valCur = interp1(x,y,posRight);
end

%%DEBUG
%hold off;
% plot(x,y);
% hold all;
% plot(posLeftMax,valLeftMax,'ro');
% plot(posLeft,valLeftMax/2,'bo')
% plot(posRightMax,valRightMax,'ro');
% plot(posRight,valRightMax/2,'bo');
%%DEBUG

width = posRight-posLeft;
xLeft = posLeft;
xRight = posRight;
yLeft =valLeftMax/2;
yRight =valRightMax/2;

% Figure for measureZ.m in the SIM analysis pipeline.
if exist('fig','var')
    figure(fig),
    subplot(2,2,2)
    plot(x,y)
    hold on, plot(x,spline)
    plot(posLeftMax,valLeftMax,'+')
    plot(posRightMax,valRightMax,'+')
    plot(xLeft,yLeft,'+')
    plot(xRight,yRight,'+')
    hold off
end

end
