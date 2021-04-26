%===============================================================%
% Daniel Hortelano Roig (11/11/2020)
% daniel.hortelanoroig@materials.ox.ac.uk 

% Checks structures in input file for inconsistencies.
%===============================================================%
%%
% Nodal network:
if any(rn(:,1) > dx) || any(rn(:,2) > dy) || any(rn(:,3) > dz)
    warning('Material volume does not fully contain the nodal network.')
end

% Drag system:
if min(glidecoefficients,[],'all') / max(glidecoefficients,[],'all') < eps(1e1)
    warning('Very large disparity of drag coefficients.')
end

% Topological constraints:
if ~(2*lmin < lmax)
    warning('Recommended constraints on lmax and lmin not imposed.')
end

if ~(0 <= areamin && 4*areamin < areamax && areamax <= lmax^2*sqrt(3)/4)
    warning('Recommended topological constraints on areamax and areamin not imposed.')
end