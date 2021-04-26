function [slipSystemsCellCart] = SlipSystemsToCartesian(slipSystemsCell, a1,a2,a3,a4)
	% INPUT:
    %       slipSystemsCell: cell of glissile slip systems -- size (numB,numPT,2)
	%		a1,a2,a3,a4: HCP lattice vectors -- each size (1,3)
	% OUTPUT:
	%		slipSystemsCellCart: cell of glissile slip systems in cartesian coordinates -- size (numB,numPT,2)
    
    [numBurgs,numPlaneTypes,~] = size(slipSystemsCell);
	slipSystemsCellCart = cell(numBurgs,numPlaneTypes,2);
    
    for b = 1:numBurgs
        for t = 1:numPlaneTypes
            
            % Burgers vector:
            burgsHCP = slipSystemsCell{b,t,1};
            burgsCart = HCPToCartesian(burgsHCP, a1,a2,a3,a4);
            slipSystemsCellCart{b,t,1} = burgsCart;
            
            % Slip plane:
            planesmb = slipSystemsCell{b,t,2};
            if isempty(planesmb)
                continue % Skip if this slip system is sessile
            end
            normalsCart = PlanesMBToCartesianNormals(planesmb, a1,a2,a3,a4);
            slipSystemsCellCart{b,t,2} = normalsCart;
        end
    end
end