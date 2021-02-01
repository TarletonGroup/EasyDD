function [climbDrag] = calcClimbDrag
% DHT
% Climb drag

global k
global MU_SI a_SI Bglide_SI
global atomVol vacFormEnergy vacMigrEnergy diffCoeffPreExpo
global radiusCore radiusInfinity
global temperature

% DDLab base units are:
% [D] = Bglide (glide drag factor = [PT] = [FT/L^2])
% [L] = a (lattice parameter)
% [P] = MU (shear modulus = [F/L^2])

% Switch to DDLab units
k_DD = k/(MU_SI*a_SI^3);
atomVol_DD = atomVol/a_SI^3;
diffCoeffPreExpo_DD = diffCoeffPreExpo*Bglide_SI/(MU_SI*a_SI^2);

% Calculate self-diffusion coefficient in DDLab units
selfDiffCoeff_DD = diffCoeffPreExpo_DD*...
                   exp(-(vacFormEnergy+vacMigrEnergy)/(k*temperature));
              
% Calculate Beclimb in DDLab units
climbDrag = k_DD*temperature*log(radiusInfinity/radiusCore)/...
            (2*pi*selfDiffCoeff_DD*atomVol_DD);

end

