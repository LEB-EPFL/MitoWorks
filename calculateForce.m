function [DaughterForce, ForceError] = calculateForce(FricCoeff,DaughterSpeed,Temp)
    
    kB = 1.38064852e-23;
% disp('Estimating pulling force...');
DaughterForce=zeros(size(DaughterSpeed,1),2);
ForceError=zeros(size(DaughterSpeed,1),2);
for j=1:2
    for i=1:size(DaughterSpeed,1)
        DaughterForce(i,j)=FricCoeff(j)*abs(DaughterSpeed(i,j)*1e-9);
        ForceError(i,j)=sqrt(2*FricCoeff(j)*kB*(273+Temp));
    end
end
end