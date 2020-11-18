%===============================================================%
% Daniel Hortelano Roig (11/11/2020)
% daniel.hortelanoroig@materials.ox.ac.uk 

% Adds all relative paths required to run the simulation.
%===============================================================%

%% Source folder paths

mainpath = '../';
addpath(mainpath); % Main
addpath(genpath([mainpath 'src/modules/'])); % Selected subprocesses
addpath(genpath([mainpath 'src/init/'])); % Preamble functions

%% Output folder paths

outputpath = append('output/',INPUTNAME,'/');
if ~exist([mainpath outputpath], 'dir')
       mkdir([mainpath outputpath]);
end