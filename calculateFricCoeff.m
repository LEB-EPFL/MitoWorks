function FricCoeff = calculateFricCoeff(Eta,D,L)


disp('Estimating friction coefficient...');
% Construct a questdlg with three options
    choice = questdlg('How would you approximate the shape of the mitochondrion?', ...
        'Shape approximation', ...
        'Sphere','Cylinder','Quit','Quit');
    % Handle response
switch choice
    case 'Sphere'
        
        FricCoeff = 6*pi()*Eta*D/2;
        
    case 'Cylinder'
        r = L/D;
        eps = (log(2*L/D))^(-1);

        if r<=1
            ShapeCorr = 1 + 0.437*r - 0.0749*r^3 + 0.0623*r^5 - 0.025*r^7;
            FricCoeff = 8*Eta*D*ShapeCorr;

        elseif r>1 && r<4
            ShapeCorr = 1.0276 + 0.3963*r - 0.0259*r^2 + 0.0014*r^3;
            FricCoeff = 8*Eta*D*ShapeCorr;
        else
            ShapeCorr = 0.0244 + 0.5504*eps + 3.328*eps^2 - 2.971*eps^3;
            FricCoeff = 2*pi()*Eta*L*ShapeCorr;
        end
    case 'Quit'
        error('Quitting...');
end

clear choice;

    
    
end