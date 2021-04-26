function [slipSystemsCellCartRot,a1Rot,a2Rot,a3Rot,a4Rot] = ...
    RotateCartesianSlipSystems(slipSystemsCellCart,a1,a2,a3,a4,rotationMatrix)
	% INPUT:
	%		slipSystemsCellCart: cartesian slip systems -- size (numBurgs,numPlaneTypes,2)
    %       a1,a2,a3,a4: HCP lattice vectors -- each size (1,3)
    %       rotationMatrix: unitary rotation matrix -- size (3,3)
	% OUTPUT:
	%		slipSystemsCellRot: rotated cartesian slip systems -- size (numBurgs,numPlaneTypes,2)
    %       a1rot,a2rot,etc: rotated HCP lattice vectors -- each size (1,3)
    
    
    %%% Rotate HCP lattic vectors
    
    a1Rot = (rotationMatrix * a1')';
    a2Rot = (rotationMatrix * a2')';
    a3Rot = (rotationMatrix * a3')';
    a4Rot = (rotationMatrix * a4')';
    
    
    %%% Rotate slip systems
    
    [numBurgs,numPlaneTypes,~] = size(slipSystemsCellCart);
	slipSystemsCellCartRot = cell(numBurgs,numPlaneTypes,2);
    
    for b = 1:numBurgs
        for t = 1:numPlaneTypes
            
            % Burgers vector:
            burgsCartOld = slipSystemsCellCart{b,t,1}; % Burgers vector
            burgsCartNew = (rotationMatrix * burgsCartOld')';
            slipSystemsCellCartRot{b,t,1} = burgsCartNew;
            
            % Slip plane:
            normalsCartOld = slipSystemsCellCart{b,t,2}; % Normal vectors
            if isempty(normalsCartOld)
                continue % Skip if this slip system is sessile
            end
            normalsCartNew = (rotationMatrix * normalsCartOld')';
            slipSystemsCellCartRot{b,t,2} = normalsCartNew;
        end
    end
end