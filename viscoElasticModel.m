function [F1, F2, T] = viscoElasticModel(Speed1proj, Speed2proj, Eta, MidPointDM1, MidPointDM2, DiametersSnakeDM1, DiametersSnakeDM2, BreakFrame,EndFrame, Pix, MinDiameter)
% calculates the tension and pulling force from the visco-elastic model by
% measuring the speed of the leading edge retraction after fission
%
%   INPUTS:
%   Speed1proj and Speed2proj = speeds of daughter mitochondria in nm
%   Eta = viscosity
%   MidPointDM1 and MidPointDM2 = backbone of daughter mitochondria
%   D1 and D2 = diameters of daughter mitochondria (from mesh) in pixels
%   DiametersSnakeDM1 and DiametersSnakeDM2 = diameters of daughter mitochondria in pixels
%   BreakFrame = last frame before fission
%   EndFrame = last frame
%   Pix = pixel size in nm
%   MinDiameter = min diameter of constriction site before fission in nm
%   
%   OUTPUTS:
%   F1 and F2 = forces exerted on daugher mitochondria 
%   T = tension at constriction site in N/nm

    % measure geometry of daughter mito
    for f=BreakFrame:EndFrame
        L1(f)=calculateBackboneLength(MidPointDM1,f)*Pix; % in nm
        L2(f)=calculateBackboneLength(MidPointDM2,f)*Pix; % in nm
        D1(f)=mean(DiametersSnakeDM1{f})*Pix; % in nm
        D2(f)=mean(DiametersSnakeDM2{f})*Pix; % in nm
    end
    L0_1=min(L1(BreakFrame:EndFrame));
    L0_2=min(L2(BreakFrame:EndFrame));
    
    % calculate force F=v*eta*pi*r^2/l_0
    F1=((Speed1proj*Eta/(L0_1)))*pi()*(D1(BreakFrame+1)*1e-9/2)^2; % in N
    F2=((Speed2proj*1e-9*Eta/(L0_2*1e-9)))*pi()*(D2(BreakFrame+1)*1e-9/2)^2; % in N

    % calculate tension T=F/(2*pi*r)
    T= abs(F1-F2)/(2*((MinDiameter(BreakFrame))*1e-9)); % in N/m

end
