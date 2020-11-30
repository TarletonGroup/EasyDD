function [mainpath] = addPaths(mainpath, varargin)
    %===============================================================%
    % Daniel Hortelano Roig (11/11/2020)
    % daniel.hortelanoroig@materials.ox.ac.uk 

    % Adds all relative paths required to run the simulation.
    %===============================================================%

    %% Source folder paths

    addpath(mainpath); % Main path
    
    for i = 1:nargin-1
        
        if ~exist([mainpath varargin{i}], 'dir')
            mkdir([mainpath varargin{i}]);
        end
        
        addpath(genpath([mainpath varargin{i}])); % Selected paths
    end
end