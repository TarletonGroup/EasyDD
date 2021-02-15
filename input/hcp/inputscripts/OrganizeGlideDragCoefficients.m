function [glideCoefficients] = OrganizeGlideDragCoefficients( ...
numBurgsA,numBurgsCA,numBurgsC,numBurgs,numPlaneTypes, ... % Organization-related inputs
glideCoeffsStruct)                                         % Glide coefficient structure
    % INPUT:
    %       numBurgsA, etc: Number of that Burgers vector type -- each size (1)
    %       numBurgs, numPlaneTypes: Number of Burgers vectors and plane types -- each size (1)
    %       glideCoeffsStruct: structure of glide coefficients
	% OUTPUT: 
	%       glideCoefficients: matrix of glide drag coefficients organized by Burgers vector, plane type, and orientation -- size (numBurgs,numPlaneTypes,2)
    
    % NOTE:
    % segments with the same orientation (screw/edge), plane type,
    % and Burgers vector type (i.e. <a>, <c+a>, or <c>)
    % will have the same glide drag coefficients.
    % glideCoefficients(P,Q,1) == Glide drag coefficient of a pure screw
    %                             segment on a slip system with Burgers 
    %                             vector index P and plane type Q.
    % glideCoefficients(P,Q,2) == Glide drag coefficient of a pure edge
    %                             segment.
    
    assert(numBurgs == numBurgsA + numBurgsCA + numBurgsC, ...
          'numBurgs != numBurgsA + numBurgsCA + numBurgsC');
        
    % Glide drag coefficients for <a> segments (Burgers vector indices 1-3):
    Ba1Screw = glideCoeffsStruct.Ba1Screw;
    Ba1Edge = glideCoeffsStruct.Ba1Edge;
    Pr1Screw = glideCoeffsStruct.Pr1Screw;
    Pr1Edge = glideCoeffsStruct.Pr1Edge;
    PyI1Screw = glideCoeffsStruct.PyI1Screw;
    PyI1Edge = glideCoeffsStruct.PyI1Edge;
    PyII1Screw = glideCoeffsStruct.PyII1Screw;
    PyII1Edge = glideCoeffsStruct.PyII1Edge;
    % Glide drag coefficients for <c+a> segments (4-9):
    Pr4Screw = glideCoeffsStruct.Pr4Screw;
    Pr4Edge = glideCoeffsStruct.Pr4Edge;
    PyI4Screw = glideCoeffsStruct.PyI4Screw;
    PyI4Edge = glideCoeffsStruct.PyI4Edge;
    PyII4Screw = glideCoeffsStruct.PyII4Screw;
    PyII4Edge = glideCoeffsStruct.PyII4Edge;
    sP4Screw = glideCoeffsStruct.sP4Screw;
    sP4Edge = glideCoeffsStruct.sP4Edge;
    % Glide drag coefficients for <c> segments (10):
    Pr10Screw = glideCoeffsStruct.Pr10Screw;
    Pr10Edge = glideCoeffsStruct.Pr10Edge;
    % Glide drag for junction segments:
    dragSessile = glideCoeffsStruct.dragsessile;
    
    glideCoefficients = zeros(numBurgs,numPlaneTypes,2);
	
    % Defining indices:
    idxA = 1;                  % 1
    idxCA = idxA + numBurgsA;  % 4
    idxC = idxCA + numBurgsCA; % 10
    
	% <a> segments
	glideCoefficients(1:idxCA-1,:,:) = repmat(reshape( ...
                               [Ba1Screw, Ba1Edge; ...
                                Pr1Screw, Pr1Edge; ...
                                PyI1Screw, PyI1Edge; ...
								PyII1Screw, PyII1Edge; ...
								dragSessile, dragSessile] ...
								,1,numPlaneTypes,2) ... % reshape: reshaping to fit into glideCoefficients
								,numBurgsA,1,1); % repmat: repeating numBurgsA times since each <a> segment has numBurgsA possible Burgers vectors
	
	% <c+a> segments
	glideCoefficients(idxCA:idxC-1,:,:) = repmat(reshape( ...
							   [dragSessile, dragSessile; ...
								Pr4Screw, Pr4Edge; ...
								PyI4Screw, PyI4Edge; ...
								PyII4Screw, PyII4Edge; ...
								sP4Screw, sP4Edge] ...
								,1,numPlaneTypes,2) ...
								,numBurgsCA,1,1);

	% <c> segments
	glideCoefficients(idxC:end,:,:) = repmat(reshape( ...
							   [dragSessile, dragSessile; ...
								Pr10Screw, Pr10Edge; ...
								dragSessile, dragSessile; ...
								dragSessile, dragSessile; ...
								dragSessile, dragSessile] ...
								,1,numPlaneTypes,2) ...
								,numBurgsC,1,1);
    
    
    % Asserts proper size of glidecoefficients:
    assert(all(size(glideCoefficients) == [numBurgs,numPlaneTypes,2]), ...
    'size(glideCoefficients) should be [numBurgs,numPlaneTypes,2], but is not!');
end