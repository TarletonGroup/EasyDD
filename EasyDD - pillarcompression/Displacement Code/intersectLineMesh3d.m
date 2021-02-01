function [points, pos, faceInds] = intersectLineMesh3d(line, vertices, faces, varargin)
%INTERSECTLINEMESH3D Intersection points of a 3D line with a mesh
%
%   INTERS = intersectLineMesh3d(LINE, VERTICES, FACES)
%   Compute the intersection points between a 3D line and a 3D mesh defined
%   by vertices and faces.
%
%   [INTERS, POS, INDS] = intersectLineMesh3d(LINE, VERTICES, FACES)
%   Also returns the position of each intersection point on the input line,
%   and the index of the intersected faces.
%   If POS > 0, the point is also on the ray corresponding to the line. 
%   
%   Example
%     [V, F] = createCube;
%     line = [.2 .3 .4 1 0 0];
%     pts = intersectLineMesh3d(line, V, F)
%     pts =
%         1.0000    0.3000    0.4000
%              0    0.3000    0.4000
%
%   See also
%   meshes3d, triangulateFaces, intersectLineTriangle3d
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2011-12-20,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.

% tolerance for detecting if a point is 
tol = 1e-12;
if ~isempty(varargin)
    tol = varargin{1};
end


% ensure the mesh has triangular faces
tri2Face = [];
if iscell(faces) || size(faces, 2) ~= 3
    [faces, tri2Face] = triangulateFaces(faces);
end

% find triangle edge vectors
t0  = vertices(faces(:,1), :);
u   = vertices(faces(:,2), :) - t0;
v   = vertices(faces(:,3), :) - t0;

% triangle normal
n   = normalizeVector3d(vectorCross3d(u, v));

% direction vector of line
dir = line(4:6);

% vector between triangle origin and line origin
w0 = bsxfun(@minus, line(1:3), t0);

a = -dot(n, w0, 2);
b = dot(n, repmat(dir, size(n, 1), 1), 2);

valid = abs(b) > tol & vectorNorm3d(n) > tol;

% compute intersection point of line with supporting plane
% If pos < 0: point before ray
% IF pos > |dir|: point after edge
pos = a ./ b;

% coordinates of intersection point
points = bsxfun(@plus, line(1:3), bsxfun(@times, pos, dir));


%% test if intersection point is inside triangle

% normalize direction vectors of triangle edges
uu  = dot(u, u, 2);
uv  = dot(u, v, 2);
vv  = dot(v, v, 2);

% coordinates of vector v in triangle basis
w   = points - t0;
wu  = dot(w, u, 2);
wv  = dot(w, v, 2);

% normalization constant
D = uv.^2 - uu .* vv;

% test first coordinate
s = (uv .* wv - vv .* wu) ./ D;
% ind1 = s < 0.0 | s > 1.0;
ind1 = s < -tol | s > (1.0 + tol);
points(ind1, :) = NaN;
pos(ind1) = NaN;

% test second coordinate, and third triangle edge
t = (uv .* wu - uu .* wv) ./ D;
% ind2 = t < 0.0 | (s + t) > 1.0;
ind2 = t < -tol | (s + t) > (1.0 + tol);
points(ind2, :) = NaN;
pos(ind2) = NaN;

% keep only interesting points
inds = ~ind1 & ~ind2 & valid;
points = points(inds, :);

pos = pos(inds);
faceInds = find(inds);

% convert to face indices of original mesh
if ~isempty(tri2Face)
    faceInds = tri2Face(faceInds);
end
function [tri, inds] = triangulateFaces(faces)
%TRIANGULATEFACES Convert face array to an array of triangular faces 
%
%   TRI = triangulateFaces(FACES)
%   Returns a 3-columns array of indices, based on the data stored in the
%   argument FACES:
%   - if FACES is a N-by-3 array, returns the same array
%   - if FACES is a N-by-4 array, returns an array with 2*N rows and 3
%       columns, splitting each square into 2 triangles (uses first and
%       third vertex of each square as diagonal).
%   - if FACES is a cell array, split each face into a set of triangles,
%       and returns the union of all triangles. Faces are assumed to be
%       convex.
%
%   [TRI INDS] = triangulateFaces(FACES)
%   Also returns original face index of each new triangular face. INDS has
%   the same number of rows as TRI, and has values between 1 and the
%   number of rows of the original FACES array.
%
%
%   Example
%     % create a basic shape
%     [n e f] = createCubeOctahedron;
%     % draw with plain faces
%     figure;
%     drawMesh(n, f);
%     % draw as a triangulation
%     tri = triangulateFaces(f);
%     figure;
%     patch('vertices', n, 'faces', tri, 'facecolor', 'r');
%
%   See also
%   meshes3d, drawMesh, mergeCoplanarFaces
%

% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2008-09-08,    using Matlab 7.4.0.287 (R2007a)
% Copyright 2008 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.

%% Tri mesh case: return original set of faces

if isnumeric(faces) && size(faces, 2) == 3
    tri = faces;
    if nargout > 1
        inds = (1:size(faces, 1))';
    end
    return;
end


%% Square faces: split each square into 2 triangles

if isnumeric(faces) && size(faces, 2) == 4
    nf = size(faces, 1);
    tri = zeros(nf * 2, 3);
    tri(1:2:end, :) = faces(:, [1 2 3]);
    tri(2:2:end, :) = faces(:, [1 3 4]);
    
    if nargout > 1
        inds = kron(1:size(faces, 1), ones(1,2))';
    end
    
    return;
end


%% Pentagonal faces (for dodecahedron...): split into 3 triangles

if isnumeric(faces) && size(faces, 2) == 5
    nf = size(faces, 1);
    tri = zeros(nf * 3, 3);
    tri(1:3:end, :) = faces(:, [1 2 3]);
    tri(2:3:end, :) = faces(:, [1 3 4]);
    tri(3:3:end, :) = faces(:, [1 4 5]);
    
    if nargout > 1
        inds = kron(1:size(faces, 1), ones(1,2))';
    end
    
    return;
end


%% Faces as cell array 

% number of faces
nf  = length(faces);

% compute total number of triangles
ni = zeros(nf, 1);
for i = 1:nf
    % as many triangles as the number of vertices minus 1
    ni(i) = length(faces{i}) - 2;
end
nt = sum(ni);

% allocate memory for triangle array
tri = zeros(nt, 3);
inds = zeros(nt, 1);

% convert faces to triangles
t = 1;
for i = 1:nf
    face = faces{i};
    nv = length(face);
    v0 = face(1);
    for j = 3:nv
        tri(t, :) = [v0 face(j-1) face(j)];
        inds(t) = i;
        t = t + 1;
    end
end
function vn = normalizeVector3d(v)
%NORMALIZEVECTOR3D Normalize a 3D vector to have norm equal to 1
%
%   V2 = normalizeVector3d(V);
%   Returns the normalization of vector V, such that ||V|| = 1. Vector V is
%   given as a row vector.
%
%   If V is a N-by-3 array, normalization is performed for each row of the
%   input array.
%
%   See also:
%   vectors3d, vectorNorm3d
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 29/11/2004.
%

% HISTORY
% 2005-11-30 correct a bug
% 2009-06-19 rename as normalizeVector3d
% 2010-11-16 use bsxfun (Thanks to Sven Holcombe)

vn   = bsxfun(@rdivide, v, sqrt(sum(v.^2, 2)));
function c = vectorCross3d(a,b)
%VECTORCROSS3D Vector cross product faster than inbuilt MATLAB cross.
%
%   C = vectorCross3d(A, B) 
%   returns the cross product of the 3D vectors A and B, that is: 
%       C = A x B
%   A and B must be N-by-3 element vectors. If either A or B is a 1-by-3
%   row vector, the result C will have the size of the other input and will
%   be the  concatenation of each row's cross product. 
%
%   Example
%     v1 = [2 0 0];
%     v2 = [0 3 0];
%     vectorCross3d(v1, v2)
%     ans =
%         0   0   6
%
%
%   Class support for inputs A,B:
%      float: double, single
%
%   See also DOT.

%   Sven Holcombe

% needed_colons = max([3, length(size(a)), length(size(b))]) - 3;
% tmp_colon = {':'};
% clnSet = tmp_colon(ones(1, needed_colons));
% 
% c = bsxfun(@times, a(:,[2 3 1],clnSet{:}), b(:,[3 1 2],clnSet{:})) - ...
%     bsxfun(@times, b(:,[2 3 1],clnSet{:}), a(:,[3 1 2],clnSet{:}));

sza = size(a);
szb = size(b);

% Initialise c to the size of a or b, whichever has more dimensions. If
% they have the same dimensions, initialise to the larger of the two
switch sign(numel(sza) - numel(szb))
    case 1
        c = zeros(sza);
    case -1
        c = zeros(szb);
    otherwise
        c = zeros(max(sza, szb));
end

c(:) =  bsxfun(@times, a(:,[2 3 1],:), b(:,[3 1 2],:)) - ...
        bsxfun(@times, b(:,[2 3 1],:), a(:,[3 1 2],:));
function n = vectorNorm3d(v)
%VECTORNORM3D Norm of a 3D vector or of set of 3D vectors
%
%   N = vectorNorm3d(V);
%   Returns the norm of vector V.
%
%   When V is a N-by-3 array, compute norm for each vector of the array.
%   Vector are given as rows. Result is then a N-by-1 array.
%
%   NOTE: compute only euclidean norm.
%
%   See Also
%   vectors3d, normalizeVector3d, vectorAngle3d, hypot3
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 21/02/2005.

%   HISTORY
%   19/06/2009 rename as vectorNorm3d

n = sqrt(sum(v.*v, 2));