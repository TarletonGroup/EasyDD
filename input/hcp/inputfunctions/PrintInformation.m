%===============================================================%
% Daniel Hortelano Roig (11/11/2020)
% daniel.hortelanoroig@materials.ox.ac.uk 

% Prints identifying information from input file.
%===============================================================%
%%
% Material
fprintf('================ %s ================\n\n', mfilename);
fprintf('Material: %s\n\n',             materialname);
fprintf('Temperature: %0.2e K\n',       temperature);
fprintf('Shear modulus: %0.2e GPa\n',   MU*pressureSI/1e9);
fprintf('Poisson ratio: %0.2e\n',       NU);
fprintf('c/a ratio: %0.2e\n',           HCPc);

% Selected subprocesses
fprintf('\n--------------   Subprocesses   --------------\n\n');
fprintf('Simulation type: %s\n',        prescribeDofs);
fprintf('Boundary condition: %s\n',     boundaryConditions);
fprintf('Save condition: %s\n',         storeBC);
fprintf('Integrator function: %s\n',    integrator);
fprintf('Mobility function: %s\n',      mobility);

% Simulation units in SI
fprintf('\n-------------- Simulation Units --------------\n\n');
fprintf('Distance: %0.2e m\n',          lengthSI);
fprintf('Time: %0.2e s\n',              timeSI);
fprintf('Velocity: %0.2e m/s\n',        velocitySI);
fprintf('Pressure: %0.2e Pa\n',         pressureSI);
fprintf('Force: %0.2e Pa m^2\n',        forceSI);
fprintf('Nodal force: %0.2e Pa m\n',    nodalforceSI);
fprintf('Drag: %0.2e Pa s\n',           dragSI);
fprintf('Temperature: %0.2e K\n',       temperatureSI);