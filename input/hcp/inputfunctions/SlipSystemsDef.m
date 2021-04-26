function [slipSystems] = SlipSystemsDef(...
    numBurgs, numPlaneTypes)
	% INPUT:
    %       numBurgs: number of predefined burgers vectors -- size (1)
	%		numPlaneTypes: number of predefined glissile plane types -- size (1)
	% OUTPUT:
	%		slipSystems: cell of glissile slip systems in MB indices -- size (numBurgs,numPlaneTypes,2)

% This script defines the available slip systems.

slipSystems = cell(numBurgs,numPlaneTypes,2);

% Type <a> Burgers vectors
% [b1](p)
slipSystems(1,1:5,1) = {1/3 * [2,-1,-1,0]};         % b1
slipSystems(1,:,2) = {[0,0,0,1]; ...  % 1Ba
					  [0,1,-1,0]; ... % 1Pr
					  [0,1,-1,1]; ... % 1PyI
					  [0,-1,1,1]; ... % 1PyII
					  []; ...         % 1sP (empty)
					  };
% [b2](p)
slipSystems(2,1:5,1) = {1/3 * [-1,-1,2,0]};         % b2
slipSystems(2,:,2) = {[0,0,0,1]; ...  % 2Ba
					  [1,-1,0,0]; ... % 2Pr
					  [1,-1,0,1]; ... % 2PyI
					  [-1,1,0,1]; ... % 2PyII
					  []; ...         % 2sP (empty)
					  };
% [b3](p)
slipSystems(3,1:5,1) = {1/3 * [-1,2,-1,0]};         % b3
slipSystems(3,:,2) = {[0,0,0,1]; ...  % 3Ba
					  [-1,0,1,0]; ... % 3Pr
					  [-1,0,1,1]; ... % 3PyI
					  [1,0,-1,1]; ... % 3PyII
					  []; ... 	      % 3sP (empty)
					  };

% Type <c+a> Burgers vectors
% [b4](p)
slipSystems(4,1:5,1) = {1/3 * [2,-1,-1,3]};         % b4
slipSystems(4,:,2) = {[]; ...         % 4Ba (empty)
					  [0,1,-1,0]; ... % 4Pr
					  [-1,0,1,1]; ... % 4PyI
					  [-1,1,0,1]; ... % 4PyII
					  [-2,1,1,2]; ... % 4sP
					  };
% [b5](p)
slipSystems(5,1:5,1) = {1/3 * [2,-1,-1,-3]};        % b5
slipSystems(5,:,2) = {[]; ... 		   % 5Ba (empty)
					  [0,1,-1,0]; ...  % 5Pr
					  [1,0,-1,1]; ...  % 5PyI
					  [1,-1,0,1]; ...  % 5PyII
					  [2,-1,-1,2]; ... % 5sP
					  };
% [b6](p)
slipSystems(6,1:5,1) = {1/3 * [-1,-1,2,3]};         % b6
slipSystems(6,:,2) = {[]; ...         % 6Ba (empty)
					  [1,-1,0,0]; ... % 6Pr
					  [1,0,-1,1]; ... % 6PyI
					  [0,1,-1,1]; ... % 6PyII
					  [1,1,-2,2]; ... % 6sP
					  };
% [b7](p)
slipSystems(7,1:5,1) = {1/3 * [-1,-1,2,-3]};        % b7
slipSystems(7,:,2) = {[]; ...         % 7Ba (empty)
					  [1,-1,0,0]; ... % 7Pr
					  [-1,0,1,1]; ... % 7PyI
					  [0,-1,1,1]; ... % 7PyII
					  [-1,-1,2,2]; ... % 7sP
					  };
% [b8](p)
slipSystems(8,1:5,1) = {1/3 * [-1,2,-1,3]};         % b8
slipSystems(8,:,2) = {[]; ...         % 8Ba (empty)
					  [-1,0,1,0]; ... % 8Pr
					  [1,-1,0,1]; ... % 8PyI
					  [0,-1,1,1]; ... % 8PyII
					  [1,-2,1,2]; ... % 8sP
					  };
% [b9](p)
slipSystems(9,1:5,1) = {1/3 * [-1,2,-1,-3]};        % b9
slipSystems(9,:,2) = {[]; ... 		   % 9Ba (empty)
					  [-1,0,1,0]; ...  % 9Pr
					  [0,1,-1,1]; ...  % 9PyI
					  [-1,1,0,1]; ...  % 9PyII
					  [-1,2,-1,2]; ... % 9sP
					  };

% Type <c> Burgers vectors
% [b10](p)
slipSystems(10,1:5,1) = {1/3 * [0,0,0,3]};          % b10
slipSystems(10,:,2) = {[]; ...         % 10Ba (empty)
					   [0,1,-1,0; ...  % 10PrI
                       -1,0,1,0; ...   % 10PrII
                        1,-1,0,0]; ... % 10PrIII
					   []; ...         % 10PyI (empty)
					   []; ...         % 10PyII (empty)
					   []; ...         % 10sP (empty)
					   };
end