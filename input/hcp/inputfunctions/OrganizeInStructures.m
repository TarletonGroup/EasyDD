%===============================================================%
% Daniel Hortelano Roig (29/11/2020)
% daniel.hortelanoroig@materials.ox.ac.uk 

% Organizes input variables in structures. These structures are
% not essential to run EasyDD; they are just for organization.
%===============================================================%
%% Initialization

scales = struct; % Stores simulation units and scales

%% Storage

scales.lengthSI = lengthSI;
scales.pressureSI = pressureSI;
scales.dragSI = dragSI;
scales.timeSI = timeSI;
scales.velocitySI = velocitySI;
scales.forceSI = forceSI;
scales.nodalforceSI = nodalforceSI;
scales.temperatureSI = temperatureSI;
scales.amag = amag;
scales.mumag = mumag;