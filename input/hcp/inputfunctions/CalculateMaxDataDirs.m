function [maxDataDirs] = CalculateMaxDataDirs(slipSystemsCell)
	% INPUT:
	%		slipSystemsCell: cell of glissile slip systems -- size (P,Q,2) for P Burgers vectors and Q plane types
	% OUTPUT:
	%		maxDataDirs: maximum number of data directions for the given slipsystems -- size (1)
	
	maxEdgeDataDirs = 0;
	
	for b = 1:size(slipSystemsCell,1)
		numEdgeDataDirs = size(cat(1,slipSystemsCell{b,:,2}),1);
		if maxEdgeDataDirs < numEdgeDataDirs
			maxEdgeDataDirs = numEdgeDataDirs;
		end
	end
	
	maxDataDirs = maxEdgeDataDirs + 1;
end