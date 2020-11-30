function [planeNormals,planeTypes] = CartesianSlipSystemsToPlaneNormalsAndTypes(slipSystems)
	% INPUT:
	%		slipSystems: cell of glissile slip systems -- size (numberBurgs,numberPlaneTypes,2)
	% OUTPUT:
	%		planeNormals: matrix of unit normal vectors of the slip planes -- size (P,4), where P is the number of slip systems in slipsystemsCell
	%		planeTypes: vector where each entry is the planeType of the corresponding plane in planeNormals -- size (P,1)
	
	[numBurgs,numPlaneTypes,~] = size(slipSystems);
	numSlipSystems = size(cat(1,slipSystems{:,:,2}),1);
	planeNormals = zeros(numSlipSystems,3);
	planeTypes = zeros(numSlipSystems,1);
	
	totSize = 0;
	% MB indices vertically concatenated in the same order as those hard-coded in slipsSystemsCell:
    for b = 1:numBurgs
        for t = 1:numPlaneTypes
            planeNormalsCur = slipSystems{b,t,2};

            if isempty(planeNormalsCur)
                continue
            end

            numGlissiles = size(planeNormalsCur,1);

            planeTypes(totSize+1:totSize+numGlissiles) = t;
            planeNormals(totSize+1:totSize+numGlissiles,:) = planeNormalsCur;

            totSize = totSize + numGlissiles;
        end
    end
    
    assert(totSize == numSlipSystems, 'totSize != numSlipSystems');
end