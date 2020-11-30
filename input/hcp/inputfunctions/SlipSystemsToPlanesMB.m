function [planesRefmb,planeTypes] = SlipSystemsToPlanesMB(slipSystemsCell)
	% INPUT:
	%		slipSystemsCell: cell of glissile slip systems -- size (numberBurgs,numberPlaneTypes,2)
	% OUTPUT:
	%		planesRefmb: matrix of MB indices -- size (P,4), where P is the number of slip systems in slipsystemsCell
	%		planeTypes: vector where each entry is the planeType of the corresponding plane in planesRefmb -- size (P,1)
	
	[numBurgs,numPlaneTypes,~] = size(slipSystemsCell);
	numSlipSystems = size(cat(1,slipSystemsCell{:,:,2}),1);
	planesRefmb = zeros(numSlipSystems,4);
	planeTypes = zeros(numSlipSystems,1);
	
	totSize = 0;
	% MB indices vertically concatenated in the same order as those hard-coded in slipsSystemsCell:
    for b = 1:numBurgs
		for t = 1:numPlaneTypes
			planesmbCur = slipSystemsCell{b,t,2};
			if isempty(planesmbCur)
				continue
			end
			numGlissiles = size(planesmbCur,1);
			
			planeTypes(totSize+1:totSize+numGlissiles) = t;
			planesRefmb(totSize+1:totSize+numGlissiles,:) = planesmbCur;
			
			totSize = totSize + numGlissiles;
		end
    end
    assert(totSize == numSlipSystems, 'totSize != numSlipSystems');
end