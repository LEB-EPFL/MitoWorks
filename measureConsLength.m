function ConsLength=measureConsLength(averageCurvaturesFWHM,StartFrame,EndFrame,ContourFWHMs1, ContourFWHMs2, N, IDX, Pix)
% calculate the length of the constriction site based on change from
% negative to positive curvature
%
%   INPUTS:
%   averageCurvaturesFWHM = average curvature per segment of mesh
%   StartFrame = first frame
%   EndFrame = frame before fission
%   ContourFWHMs1 and 2 = FWHM contours [x y] of both sides
%   N = size of moving average window (N/2 on each side)
%   IDX = position of FWHM min diameter
%   Pix = size of pixel in nm
%
%   OUTPUTS:
%   ConsLength = length of constriction site in 


BB=cell(length(StartFrame:EndFrame),1);
MAcurv=cell(length(StartFrame:EndFrame),1);
ConsLength=zeros(1,length(StartFrame:EndFrame));
for f=StartFrame:EndFrame
    AvgCurv=averageCurvaturesFWHM{f};
    AvgCurv(AvgCurv<0)=-1;
    AvgCurv(AvgCurv>0)=1;
    BB{f}=[(ContourFWHMs1{f}(:,1)+ContourFWHMs2{f}(:,1))/2 (ContourFWHMs1{f}(:,2)+ContourFWHMs2{f}(:,2))/2];
    AvgCurvOld=zeros(size(AvgCurv));
    while AvgCurvOld~=AvgCurv
        AvgCurvOld=AvgCurv;
    for i=1:length(averageCurvaturesFWHM{f})
        if i<=(N/2)
            MAcurv{f}(i)=mean(averageCurvaturesFWHM{f}(1:round(i+(N/2)),1));
            AvgCurv(i)=mean(AvgCurv(1:round(i+(N/2))));
        elseif i>=(length(averageCurvaturesFWHM{f})-(N/2))
            MAcurv{f}(i)=mean(averageCurvaturesFWHM{f}(round(i-(N/2):end),1));
            AvgCurv(i)=mean(AvgCurv((round(i-(N/2):end))));
        else
            MAcurv{f}(i)=mean(averageCurvaturesFWHM{f}(round(i-(N/2)):round(i+(N/2)),1));
            AvgCurv(i)=mean(AvgCurv(round(i-(N/2)):round(i+(N/2))));
        end
    end
    end
    L=0; R=0; Li=IDX(f); Ri=IDX(f);
    while L==0
        Li=Li-1;
        if Li>0
            if AvgCurv(Li)>0
                L=1;
            end
        else
            Li=1;
            L=1;
        end
            
    end
    while R==0
        Ri=Ri+1;
        if Ri<=length(AvgCurv)
            if AvgCurv(Ri)>0
                R=1;
            end
        else
            Ri=length(BB{f});
            R=1;
        end
    end
    ConsLength(f)=0;
    for j=Li:Ri-1
        ConsLength(f)=ConsLength(f)+(sqrt((BB{f}(j+1,1)-BB{f}(j,1))^2 + (BB{f}(j+1,2)-BB{f}(j,2))^2))*Pix;
    end
end


end
