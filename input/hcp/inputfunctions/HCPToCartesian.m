function [refCart] = HCPToCartesian(refmb,a1,a2,a3,a4)
	% INPUT:
	%		refmb: matrix of vectors expressed in terms of four HCP lattice vectors -- size (Q,4) for Q vectors
	% 		a1,a2,a3,a4: HCP lattice vectors -- each size (1,3)
	% OUTPUT: 
	%		refCart: matrix of corresponding vectors in cartesian basis -- size (Q,3)
	
    refCart = zeros(size(refmb,1),3);
    
	for i = 1:size(refmb,1)
	    refCart(i,:) = refmb(i,1)*a1 + refmb(i,2)*a2 + refmb(i,3)*a3 + refmb(i,4)*a4;
	end
end