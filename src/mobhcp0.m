%===============================================================%
% Daniel Hortelano Roig (30/11/2020)
% daniel.hortelanoroig@materials.ox.ac.uk

% Linear HCP mobility law function.
% Solves for v:
% f = B*v

% f: nodal force
% B: drag matrix
% v: nodal velocity

% For an overview of the HCP mobility law and notation used in this script,
% see Aubry et al: http://dx.doi.org/10.1016/j.jmps.2016.04.019
%===============================================================%

function [vn,fn] = mobhcp0(fseg,rn,links,connectivity,nodelist,conlist, ...
mobstruct, varargin)
% OUTPUT:
%		fn: matrix of nodal forces -- size (N,3) for N input nodes
%       vn: matrix of corresponding nodal response velocities -- size (N,3)

%%% Check whether number of inputs is valid:
if size(varargin,2) > 0
    % This function is not designed for this. Assume first extra input
    % argument is a rotation matrix.
    
    rotMatrix = varargin{1};
    
    rotateCoords = false;
    if ~isempty(rotMatrix)
        rotateCoords = true;
        rn(:, 1:3) = rn(:, 1:3) * rotMatrix;
        fseg(:, 1:3) = fseg(:, 1:3) * rotMatrix;
        fseg(:, 4:6) = fseg(:, 4:6) * rotMatrix;
        links(:, 3:5) = links(:, 3:5) * rotMatrix;
        links(:, 6:8) = links(:, 6:8) * rotMatrix;
    end
end

%% Mobility structure extraction %%

%%% Glissile cartesian slip systems:
slipsystems = mobstruct.slipsystemsref; % Cartesian
[numburgs,numplanetypes,~] = size(slipsystems);

%%% Reference Burgers vectors:
burgscart = mobstruct.burgscart; % Cartesian

%%% Drag coefficients:
dragline = mobstruct.dragline; % Line
dragclimb = mobstruct.dragclimb; % Climb
glidecoefficients = mobstruct.glidecoefficients; % Glide
dragsessile = mobstruct.dragsessile; % Sessile ("junction")

%%% Critical resolved shear stress (CRSS):
crsscoefficients = mobstruct.crsscoefficients;
crsssessile = mobstruct.crsssessile;

%%% Dimension checks:
assert(all(size(burgscart) == [numburgs,3]), ...
    'Assertion failed: all(size(burgscart) == [numburgs,3])')
assert(all(size(glidecoefficients) == [numburgs,numplanetypes,2]), ...
    'Assertion failed: all(size(glidecoefficients) == [numburgs,numplanetypes,2])')
assert(all(size(crsscoefficients) == [numburgs,numplanetypes + 1]), ...
    'Assertion failed: all(size(crsscoefficients) == [numburgs,numplanetypes + 1])')


%% Mobility law parameters %%

maxe = 2; % Max number of pure edge segment decompositions

%%% Struct containing numerical tolerances for mobility law:
mobeps = struct(...
    'eps_seglen', 1e7, ... % Max segment length before linprog could have problems
    'eps_reflist', 1e-2, ... % Min difference between actual and reference vectors
    'eps_decomp', 1e-4, ... % Min relative magnitude of decomposed segment component
    'eps_rconddrag', 1e-9, ... % Min drag matrix inverse condition number
    'eps_dragline', 1/dragline, ... % Drag line coefficient to improve conditioning
    'eps_noise', 1e-5, ... % Noise parameter to improve conditioning
    'eps_noisemax', 1e-2 ... % Max allowed value for eps_rcondnoise
    );


%% Mobility law preprocessing %%

%%% Identify the list and connectivity of the nodes for which the velocity will be computed:
if isempty(nodelist) || isempty(conlist)
    [nodelist,conlist] = SetupNodelistAndConlist(rn,connectivity);
    % NOTE:
    % conlist(i,:) = [M 1 2 ... (M-1) M 0 ... 0 ], where M = number of
    %                segments connected to node i
    % nodelist = [1 2 3 ... size(rn,1)]
end
numnodes = size(nodelist,1);

%%% Preallocate nodal force and velocity:
fn = zeros(numnodes,3);
vn = zeros(numnodes,3);


%% MOBILITY LAW %%


for n = 1:numnodes
    %% Determine dislocation properties %%
    
	n0 = nodelist(n); % nth node in nodelist
	numNbrs = conlist(n,1); % Number of segments connected to node n0
	fn(n,:) = zeros(1,3);
    
	burg = zeros(numNbrs,3); % Burgers vectors
	ndir = zeros(numNbrs,3); % Normal vectors
	linedir = zeros(numNbrs,3); % Line vectors
	
	seglenvec = zeros(numNbrs,3); % Length vectors
	seglen = zeros(numNbrs,1); % Length scalars
	
	% Segment length decomposition information
	screwseglen = zeros(numNbrs,1); % Screw length component
    unitburg = zeros(numNbrs,3); % Screw unit vector
	edgeseglen = zeros(numNbrs,maxe); % Edge length components
	edgeseglenunitvec = zeros(numNbrs,maxe,3); % Edge unit vectors
	
	% Glide drag coefficients
	screwdragcoeff = zeros(numNbrs,1);
	edgedragcoeff = zeros(numNbrs,maxe);
	
	junctionexists = zeros(numNbrs,1); % =1: junction, =0: otherwise
    
    crsspassed = zeros(numNbrs,1); % =0: segment force overcame CRSS, =1: otherwise
    
	% Loop over all segment i connected to n0:
	for i = 1:numNbrs
		%% Determine connected segment properties %%
        
        % ith segment connected to n0 in conlist:
		ii = conlist(n,i+1);
        % segment-id of segment i:
		linkid = connectivity(n0,2*ii); 
        % "position" of n0 in segment i (= 1 or = 2):
		posinlink = connectivity(n0,2*ii+1);
        % node-id of the other node in segment i:
		n1 = links(linkid,3-posinlink);
        
		% Calculate the length of segment i:
		seglenvec(i,:) = rn(n1,1:3) - rn(n0,1:3);
		seglen(i) = norm(seglenvec(i,:));
		
        % Ignore segments with zero length:
        if seglen(i) == 0.0
			continue
        end
        
		% Dislocation properties of segment i:
		burg(i,:) = links(connectivity(n0,2*ii),3:5);
		ndir(i,:) = links(connectivity(n0,2*ii),6:8); % Predefined
		linedir(i,:) = seglenvec(i,:) ./ seglen(i);
        
        
        %% Determine segment type %%
        
        % Default parameters for sessile segments:
        planetype = 6;
        pindex = -1;
        
        % Match the Burgers vector to one in the reference list:
        bindex = DetermineIndexOfVectorInReferenceList(burg(i,:), burgscart, mobeps);
        % If at least one match was found, bindex is a row vector
        % containing all matched indices as integers between 1 and
        % size(burgscart,1). If no match was found, bindex = [].
        if size(bindex,2) > 1
            bindex = bindex(1); % Select the first matching index
            warning('The Burgers vector matched more than one vector in the reference list.');
        elseif isempty(bindex)
            bindex = -1;
        end
        
        if (bindex >= 1 && bindex <= numburgs)
            % The Burgers vector is in the reference list.
			
            % All possible glissile plane MB indices and their
            % corresponding planetypes for the current Burgers vector:
            [glissilenormals,glissileplanetypes] = CartesianSlipSystemsToPlaneNormalsAndTypes(slipsystems(bindex,:,:));
            
            numglissilenormals = size(glissilenormals,1);
            
            % Match the normal vector to one in the reference list:
            pindex = DetermineIndexOfVectorInReferenceList(ndir(i,:), glissilenormals, mobeps);
            if size(pindex,2) > 1
                pindex = pindex(1); % Select the first matching index
                warning('The normal vector matched more than one vector in the reference list.');
            elseif isempty(pindex)
                pindex = -1;
            end
            
            if (pindex >= 1 && pindex <= numglissilenormals)
                % The normal vector is in the reference list.
                
                % Glissile planetype (integer between 1 and 5):
                planetype = glissileplanetypes(pindex);
                % Reassign normal vector:
                ndir(i,:) = glissilenormals(pindex,:);
                
                if planetype == 6
                    warning('Contradiction: match found in reference list for segment normal vector, but the planetype is 6.');
                end
            end
        end
		
        
        %% Determine segment decomposition and glide drag coefficients %%
        
        if planetype == 6
            % The segment is a junction.
            
            if bindex ~= -1 && pindex ~= -1
                fprintf('bindex: '); disp(bindex);
                fprintf('\npindex: '); disp(pindex);
                warning('Sessile segment; neither bindex nor pindex are -1, but at least one should be.');
            end
            
			junctionexists(i) = 1;
        else
            % The segment is NOT a junction; it is glissile.
            % NOTE: junctionexists(i) = 0 by default.
            
            if abs(seglen(i)) > mobeps.eps_seglen
                warning('Very large segment length of: %e.\n',seglen(i));
            end
            
            % Decompose segment:
            [screwseglen(i), unitburg(i,:), edgeseglen(i,:), edgeseglenunitvec(i,:,:), ...
                edgeidx1stmax, edgeidx2ndmax] = DecomposeSegment(...
                glissilenormals, pindex, burgscart(bindex,:), seglenvec(i,:));
			
            % Determine screw and edge glide coefficients:
            [screwdragcoeff(i), edgedragcoeff(i,:)] = DetermineGlideCoefficients(...
                glidecoefficients, glissileplanetypes, ...
                bindex, planetype, edgeidx1stmax, edgeidx2ndmax);
        end
        
        
        %% Determine nodal forces %%
        
        % Contribution to total nodal force on node n0 from segment i:
		fsegn0 = fseg(linkid, 3 * (posinlink - 1) + (1:3));
        
        % Superpose fsegn0 to the total force on n0:
        fn(n,:) = fn(n,:) + fsegn0;
        
        % Contribution to total nodal force on node n1 from segment i:
        posinlink1 = (posinlink == 1) + 1;
        fsegn1 = fseg(linkid, 3 * (posinlink1 - 1) + (1:3));
        if ~((posinlink == 1) || (posinlink == 2))
            error('posinlink should be 1 or 2, but is %d', ...
                posinlink);
        end
        
        fsegmid = 0.5 * (fsegn0 + fsegn1); % Nodal force on midpoint of segment i
        
        % Determine whether segment passes CRSS:
        [crsspassed(i), junctionexists(i)] = DetermineCRSSActivation(...
            crsscoefficients, crsssessile, ...
            bindex, planetype, ...
            fsegmid, ndir(i,:), seglen(i));
        % NOTE: crsspassed has no functional purpose; it's for debugging.
        
        
	end % Repeat for all connected segments
    
    
    %% Construct total drag matrix %%
    
    [Btotal] = BuildTotalDragMatrix(...
        numNbrs, junctionexists, ...
        seglen, screwseglen, edgeseglen, ...
        linedir, unitburg, ndir, edgeseglenunitvec, ...
        edgedragcoeff, screwdragcoeff, dragline, dragsessile, dragclimb, ...
        mobeps);
    
    %%% Verify properties of the drag matrix:
    if ~isreal(Btotal)
        warning('The drag matrix should be real, but is not.');
    end
    if ~issymmetric(Btotal)
        warning('The drag matrix should be symmetric, but is not.');
    end
    assert(all(size(Btotal) == [3 3]), ...
        'Assertion failed: all(size(Btotal) == [3 3])')
    
    
    %% Analyse and solve drag system
        
    %%% Linear solution:
    [vL] = LinearSystemSolver(...
        numNbrs, ...
        Btotal, fn(n,:), linedir, ...
        mobeps);
    
    %%%%%%%%%%%%%%
    vn(n,:) = vL; % Linear system
	%%%%%%%%%%%%%%
    
    if size(varargin,2) > 0
        if rotateCoords
            vn = vn * rotMatrix';
            fn = fn * rotMatrix';
        end
    end
    
end % Repeat for each node in nodelist
end % Mobility function complete



%% FUNCTION DEFINITIONS %%


%% Preprocess list of nodes  %%

function [nodeList,conList] = SetupNodelistAndConlist(rn,connectivity)

	L1 = size(rn,1); % Length of the nodelist for which the velocity will be calculated
	nodeList = linspace(1,L1,L1)';
	[L2,L3] = size(connectivity); % L2 = number of nodes; L3 = (2*maxconnectivity + 1)
    
    assert(L1 == L2, 'size(rn,1) should be equal to size(connectivity,1), but is not.');
	
	conList = zeros(L2,(L3-1)/2+1); % zeros(L2, maxconnectivity+1)
	conList(:,1) = connectivity(:,1);
	for i = 1:L2
		connumb = conList(i,1);
		conList(i,2:connumb+1) = linspace(1,connumb,connumb);
	end
end


%% Match vectors to reference list  %%

function [vIndex] = DetermineIndexOfVectorInReferenceList(...
    cartVec, cartRefVecs, mobEps)
	% INPUT:
	%		cartVec: input vector -- size (1,3)
	% 		cartRefVecs: matrix of vectors in reference list -- size (P,3)
    %       mobEps: numerical tolerances -- struct
	% OUTPUT:
	%		vIndex: vector of indexes of all matched vectors in reference list -- size (1,Q)
    
    % Determines collinearity match via cross product.
    % NOTE: if no matches are found, vIndex = [].
    
    eps_RefList = mobEps.eps_reflist; % Scalar tolerance for whether cartVec is in reference list
        
    numRefVecs = size(cartRefVecs,1);
    unitCartVec = cartVec/norm(cartVec);
    
    vIndex = zeros(1,numRefVecs); % Preallocate for max possible matches
    for v = 1:numRefVecs
        unitCartRefVec = cartRefVecs(v,:)/norm(cartRefVecs(v,:));
        if norm(cross(unitCartVec,unitCartRefVec)) < eps_RefList
            
            % If a match is found within tolerance, vIndex(v) is the index
            % (between 1 and numRefVecs) of the corresponding matching
            % value:
            vIndex(v) = v;
        end
    end
    
    % Remove zero entries from vIndex:
    vIndex = vIndex(vIndex ~= 0);
    
    if ~(isempty(vIndex) || (all(vIndex >= 1) && all(vIndex <= numRefVecs)))
        warning('vIndex contains at least one index out of the bounds of the reference list.');
    end
end


%% Transform slip systems %%

function [planeNormals, planeTypes] = CartesianSlipSystemsToPlaneNormalsAndTypes(...
    slipSystems)
	% INPUT:
	%		slipSystems: cell of glissile slip systems -- size (numberBurgs,numberPlaneTypes,2)
	% OUTPUT:
	%		planeNormals: matrix of unit normal vectors of the slip planes -- size (P,4), where P is the number of slip systems in slipsystemsCell
	%		planeTypes: vector where each entry is the planeType of the corresponding plane in planeNormals -- size (P,1)
	
	[numBurgs,numPlaneTypes,~] = size(slipSystems);
	numSlipSystems = size(cat(1,slipSystems{:,:,2}),1);
	planeNormals = zeros(numSlipSystems,3);
	planeTypes = zeros(numSlipSystems,1);
	
	totSize = 0; % Initial integer for size counter
	% MB indices vertically concatenated in the same order as those hard-coded in slipSystems:
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


%% Decompose dislocation line segment %%

function [screwSegLen, screwSegLenUnitVec, edgeSegLens, edgeSegLensUnitVec, ...
    edgeProjIdx1stMax, edgeProjIdx2ndMax] = DecomposeSegment(...
    glissileNormals, pIndex, burgVec, segLenVec)
	% INPUT:
	%		glissilenormals: glissile plane normal vectors for the current segment -- size (P,3) for P glissile planes
	% 		burgVec: Burgers vector of the segment -- size (1,3)
	% 		segLenVec: length vector of the segment -- size (1,3)
    % 		pIndex: segment index in glissileNormals  -- size (1)
	% OUTPUT:
	%		screwSegLen: screw component length -- size (1)
    %		screwSegLenUnitVec: screw component unit vector -- size (1,3)
	%		edgeSegLens: decomposed edge lengths -- size (1,2)
	%		edgeSegLensUnitVec: decomposed edge unit vectors -- size (1,2,3)
    %		edgeIdx1stMax: index corresponding to 1st largest edge projection magnitude -- size (1)
    %		edgeProjIdx2ndMax: index corresponding to 2nd largest edge projection magnitude -- size (1)
    
    % The segment length vector is decomposed into a pure screw and a pure edge
    % component. The pure edge component is further decomposed into two
    % glissile edge line directions which minimise its L1-norm.
    
    % Pure screw vector:
	burgUnitVec = burgVec/norm(burgVec);
    
    % Pure edge vector:
    edgeSegVec = segLenVec - dot(segLenVec, burgUnitVec) * burgUnitVec;
    
    %%% Find 1st largest edge projection magnitude
    
    glissileIdxs = 1:size(glissileNormals,1);
    
    % This search can be avoided by assuming edgeProjIdx1stMax = pIndex,
    % but is useful for debugging purposes.
    [edgeProjIdx1stMax, edgeProj1stMax, glissileEdgeUnitVec1stMax] = GetLargestEdgeSegVecProjection(...
        glissileNormals, glissileIdxs, burgUnitVec, edgeSegVec);
    
    %%% Find 2nd largest edge projection magnitude
        
    % Remove index of 1st max edge projection:
    glissileIdxs2ndMax = glissileIdxs(glissileIdxs ~= edgeProjIdx1stMax);
    
    % Find 2nd max edge projection:
    [edgeProjIdx2ndMax, edgeProj2ndMax, glissileEdgeUnitVec2ndMax] = GetLargestEdgeSegVecProjection(...
        glissileNormals, glissileIdxs2ndMax, burgUnitVec, edgeSegVec);
    
    %%% Checks
    
    if edgeProjIdx1stMax ~= pIndex
        warning('edgeProjIdx1stMax, %d, is NOT equal to pIndex, %d. This indicates a possible issue with the input slip systems or segment normal vector.', ...
            edgeProjIdx1stMax,pIndex);
    end
    if abs(edgeProj1stMax) == abs(edgeProj2ndMax)
        warning('edgeProj1stMax, %e, is equal to edgeProj2ndMax, %e.', ...
            abs(edgeProj1stMax),abs(edgeProj2ndMax));
    end
    
    %%% Express edge segment vector in a new basis
    
    % Express edgeSegVec in a basis composed of the 1st and 2nd largest
    % edge projections on the glissile edge line directions:
    [edgeSegLens, ~] = DecomposeVectorInBasis(...
        edgeSegVec, glissileEdgeUnitVec1stMax, glissileEdgeUnitVec2ndMax);
    
    edgeSegLens = edgeSegLens(1:2); % Select first two entries
    
    %%% Prepare output
    
    screwSegLen = dot(segLenVec,burgUnitVec);
    screwSegLenUnitVec = reshape(burgUnitVec, 1,3);
    edgeSegLens = reshape(edgeSegLens, 1,2);
    edgeSegLensUnitVec = [glissileEdgeUnitVec1stMax; ...
                         glissileEdgeUnitVec2ndMax];
    edgeSegLensUnitVec = reshape(edgeSegLensUnitVec, 1,2,3);
end


function [edgeProjIdxMax, edgeProjMax, glissileEdgeUnitVecMax] = GetLargestEdgeSegVecProjection(...
    glissileNormals, glissileIdxs, unitBurgVec, edgeSegVec)
    
    % NOTE: plane indices in glissileIdxs control which glissileNormals
    % will be compared.
    
    % Default initial guesses:
    edgeProjMax = 0;
    edgeProjIdxMax = 1;
    glissileEdgeVecMax = cross(unitBurgVec,glissileNormals(edgeProjIdxMax,:));
    glissileEdgeUnitVecMax = glissileEdgeVecMax / norm(glissileEdgeVecMax);
    
	for p = glissileIdxs
        
		normal = glissileNormals(p,:);
		glissileEdgeVec = cross(unitBurgVec,normal);
		glissileEdgeUnitVec = glissileEdgeVec/norm(glissileEdgeVec);
        
        edgeProj = dot(glissileEdgeUnitVec, edgeSegVec);
        
        if abs(edgeProj) > abs(edgeProjMax)
            edgeProjIdxMax = p;
            edgeProjMax = edgeProj;
            glissileEdgeUnitVecMax = glissileEdgeUnitVec;
        end
	end
end


function [outVec, outDualVec] = DecomposeVectorInBasis(...
    inVec, varargin)
    % INPUT:
	%		inVec: vector expressed in standard basis -- size (1,3)
	%       a1,a2,a3: three linearly independent basis vectors -- each size (1,3)
	% OUTPUT:
	%		outVec: inVec expressed in (a1,a2,a3) basis -- size (1,3)
    %		outDualVec: inVec expressed in (a1,a2,a3) dual basis -- size (1,3)
    
    % Obtain basis vectors from variable number of function inputs:
    a1 = varargin{1};
    a2 = varargin{2};
    if size(varargin,2) == 3
        a3 = varargin{3}; % Third basis vector given as input
    elseif size(varargin,2) == 2
        a3 = cross(a1,a2); % Third basis vector generated
    else
        warning('Number of basis vectors given: %d. Needed: 2 or 3.', ...
            size(varargin,2));
    end
    directBasis = [a1;a2;a3];
    
	volLat = det(directBasis); % Volume generated by basis vectors
	
	% Checking linear independence of basis vectors:
    eps_VolLat = 1e-12;
	if abs(volLat) < eps_VolLat
        warning('The volume generated by the three basis vectors is %e, which is less than the tolerance %e.\n', ...
            abs(volLat),eps_VolLat)
	end
	
	% Dual basis vectors:
	g1 = cross(a2,a3) ./ volLat;
	g2 = cross(a3,a1) ./ volLat;
	g3 = cross(a1,a2) ./ volLat;
    dualBasis = [g1;g2;g3];
    
    % NOTE:     outVec(i) = dot(g(i),inVec);
    %           outDualVec(i) = dot(a(i),inVec)
    
    % Direct basis decomposition:
    outVec = (dualBasis * inVec')';
    
    % Dual basis decomposition:
    outDualVec = (directBasis * inVec')';
end


%% Determine Segment Glide Drag Coefficients and CRSS Activation

function [screwDragCoeff, edgeDragCoeff] = DetermineGlideCoefficients(...
    glideCoefficients, glissilePlaneTypes, ...
    bIndex, planeType, edgeIdx1stMax, edgeIdx2ndMax)
    
    %%% Screw glide drag coefficients
    
    screwDragCoeff = glideCoefficients(bIndex,planeType,1);

    %%% Edge glide drag coefficients
    
    % Edge plane corresponding to 1st largest edge projection:
    edgePlaneType1stMax = glissilePlaneTypes(edgeIdx1stMax);
    edgeDragCoeff(1) = glideCoefficients(bIndex,edgePlaneType1stMax,2);
    
    % Edge plane corresponding to 2nd largest edge projection:
    edgePlaneType2ndMax = glissilePlaneTypes(edgeIdx2ndMax);
	edgeDragCoeff(2) = glideCoefficients(bIndex,edgePlaneType2ndMax,2);
end


function [crssPassed, junctionExists] = DetermineCRSSActivation(...
    crssCoefficients, crssSessile, ...
    bIndex, planeType, ...
    fSegMid, nDir, segLen)
    
    % Determines whether the nodal force (units: N/m) on a segment is
    % sufficient to overcome the CRSS. First, the CRSS of the segment is
    % determined based on the type of plane it is gliding on. Next, the
    % segment force is compared to the CRSS. If the CRSS is NOT overcome,
    % crssPassed = 0 and the segment is treated as a junction (i.e.
    % sessile) and junctionExists = 1.
    
    % Determine CRSS:
    if planeType == 6
        crss = crssSessile; % Sessile CRSS
    else
        crss = crssCoefficients(bIndex,planeType);
    end
    
    % Impose CRSS:
    fSegn0Plane = fSegMid - dot(fSegMid,nDir) * nDir;
    if norm(fSegn0Plane) / segLen < crss
        % Segment force was not sufficient to overcome CRSS.
        crssPassed = 0;
        % Treat this segment as a junction:
        junctionExists = 1;
    else
        % Segment force was sufficient to overcome CRSS.
        crssPassed = 1;
        junctionExists = 0;
    end
end


%% Build drag matrix %%

function [totalDragMat] = BuildTotalDragMatrix(...
    numConnectionsToNode, junctionExists, ...
    segLen, screwSegLen, edgeSegLen, ...
    lineDir, unitBurg, nDir, edgeSegLenUnitVec, ...
    edgeDragCoeff, screwDragCoeff, dragLine, dragSessile, dragClimb, ...
    mobEps)
	
    % Extraction:
    eps_Decomp = mobEps.eps_decomp;
    
	totalDragMat = zeros(3,3); % Initialisation
    
    for i = 1:numConnectionsToNode
        
        % Segments with zero length do not contribute:
        if segLen(i) == 0
            continue
        end
        
        % Junction contribution:
        if junctionExists(i) == 1
            
            % Construct the junction drag matrix:
            junctDragMatrix = BuildJunctionDragMatrix(...
                lineDir(i,:), dragLine, dragSessile);
            
            % Add contribution to total drag function:
            totalDragMat = totalDragMat + ...
                CalculateDragMatrixContribution(segLen(i), junctDragMatrix);
            
            continue
        end

        % Screw contribution:
        if abs(screwSegLen(i) / segLen(i)) > eps_Decomp
            
            % Unit direction vectors:
            lineVec = unitBurg(i,:);
            climbVec = nDir(i,:);
            glideVec = cross(lineVec,nDir(i,:))/norm(cross(lineVec,nDir(i,:)));
            
            % Construct the screw drag matrix:
            screwDragMat = BuildGlissileDragMatrix(...
                lineVec, climbVec, glideVec, dragLine, dragClimb, screwDragCoeff(i));
            
            totalDragMat = totalDragMat + ...
                CalculateDragMatrixContribution(screwSegLen(i), screwDragMat);
        end
        
        % Edge contribution:
        if abs(norm(edgeSegLen(i,:)) / segLen(i)) > eps_Decomp
            
            for p = 1:2
                
                % Unit direction vectors:
                lineVec = reshape(edgeSegLenUnitVec(i,p,:),1,3);
                climbVec = cross(unitBurg(i,:),lineVec)/norm(cross(unitBurg(i,:),lineVec));
                glideVec = unitBurg(i,:);
                
                % Construct the edge drag matrix:
                edgeDragMat = BuildGlissileDragMatrix(...
                    lineVec, climbVec, glideVec, dragLine, dragClimb, edgeDragCoeff(i,p));
                
                totalDragMat = totalDragMat + ...
                    CalculateDragMatrixContribution(edgeSegLen(i,p), edgeDragMat);
            end
        end
    end
end


function [dragMatrixContribution] = CalculateDragMatrixContribution(segLength, dragMatrix)
	
	dragMatrixContribution = 0.5 * abs(segLength) .* dragMatrix;
end


function [junctDragMatrix] = BuildJunctionDragMatrix(lineVec, junctLineDragCoeff, junctGlideDragCoeff)
	
    sysDim = size(lineVec,2); % = 3 for full linear system
    
	junctLineDirMat  = lineVec'  * lineVec;
	junctGlideDirMat = eye(sysDim) - junctLineDirMat;
	
	junctDragMatrix = junctLineDragCoeff * junctLineDirMat + junctGlideDragCoeff * junctGlideDirMat;
end


function [glissileDragMat] = BuildGlissileDragMatrix(...
    lineVec, climbVec, glideVec, lineDragCoeff, climbDragCoeff, glideDragCoeff)

	lineDirMat = lineVec'  * lineVec;
	climbDirMat = climbVec' * climbVec;
	glideDirMat = glideVec' * glideVec;
	
	glissileDragMat = lineDragCoeff * lineDirMat + ...
        climbDragCoeff * climbDirMat + ...
        glideDragCoeff * glideDirMat;
end


%% Solve and analyse linear system %%

function [vL] = LinearSystemSolver(...
    numConnectionsToNode, ...
    Btot, fVec, lineDir, ...
    mobEps)
    % INPUT:
	%		Btot: total drag matrix -- size (3,3)
	% 		fVec: nodal force vector -- size (1,3)
    %		lineDir: line directions of Q connected links -- size (Q,3)
    % 		mobEps: numerical tolerances -- struct
	% OUTPUT:
	%		vL: solution to linear approximation -- size (1,3)
    
    % Extraction:
    eps_DragLine = mobEps.eps_dragline;
    eps_rCondDrag = mobEps.eps_rconddrag;
    
    fVec = fVec'; % Changed to size (3,1)
    
    %%% Prepare for diagnosis and potential restructuring of linear system:
    % Diagnostic information:
    [eVecs,eVals] = eigs(Btot); % Eigenvectors/eigenvalues
    imageIdxs = (1:size(eVals,1))'; % (= [1 2 3]') % Image space indices
    % Restructuring matrices:
    BtotL = Btot;
    fVecL = fVec;
    
    %%% Analyse linear system singularities:
    if any(diag(eVals) < eps(10))
        % Drag matrix contains at least one zero eigenvalue to machine
        % precision; corresponding eigenvectors are nullspace basis
        % vectors. The linear system will be solved in its eigenbasis
        % with reduced dimensionality (nullspace dimensions removed) to
        % improve stability.
        
        % NOTE: this disables nodal motion in nullspace direction vectors.
        
        % Set near-zero eigenvalues to zero:
        eVals(abs(eVals) < eps(10)) = 0; % eps(10) \approx 1.78e-15
        
        imageIdxs = (diag(eVals) ~= 0); % Nonzero eigenvalue indices
        
        % Rotate nodal force vector into eigenbasis:
        fVecL = eVecs * fVecL;
        
        % Remove dimensions corresponding to zero eigenvalues:
        fVecL = fVecL(imageIdxs);
        BtotL = eVals(imageIdxs,imageIdxs);
    end
    dimSys = size(Btot,1); % Original system dimension
    dimSysL = size(BtotL,1); % Restructured system dimension
    
    %%% Analyse drag matrix conditioning:
    if rcond(BtotL) < eps_rCondDrag
        % Drag matrix is poorly conditioned; try decreasing the disparity
        % of the magnitude of its set of eigenvalues.
        
        junctGlideDragCoeff = 0;
        % Increase drag in the line directions of the connected segments:
        for i = 1:numConnectionsToNode
            
            segLineVec = lineDir(i,:)';
            if dimSysL < dimSys
                
                % Rotate line vector into eigenbasis:
                segLineVec = eVecs * segLineVec;
                % Reduce dimension:
                segLineVec = segLineVec(imageIdxs);
            end
            BtotL = BtotL + ...
                BuildJunctionDragMatrix(segLineVec', eps_DragLine, junctGlideDragCoeff);
        end
        
        if rcond(BtotL) < eps_rCondDrag
            % Linear system is still poorly conditioned; increase the
            % stiffness of the node in isotropic directions by adding
            % unform isotropic noise (i.e. scaled identity matrix).
            
            BtotL = BtotL + BuildNoiseMatrix(BtotL, mobEps);
        end
    end
    
    %%% Final nodal velocity
    
    vL = BtotL \ fVecL;
    
    if dimSysL < dimSys
        
        % Replace the nullspace dimension with zeros:
        vL0 = zeros(3,1);
        vL0(imageIdxs) = vL;
        % Undo the rotation into the eigenbasis:
        vL0 = eVecs' * vL0; % Note orthogonality of eVecs
        
        vL = vL0;
    end
    
    % Check consistency of system size:
    assert(size(Btot,1) == size(vL,1), ...
        'Assertion failed: size(Btot,1) == size(vnL,1)')
    assert(size(BtotL,1) == size(fVecL,1), ...
            'Assertion failed: size(BtotL,1) == size(fnL,1)')
    assert(size(imageIdxs,1) == size(fVecL,1), ...
            'Assertion failed: size(imageIdxs,1) == size(fnL,1)')
    
    vL = vL'; % Changed to size (1,3)
    
    % Check validity of nodal velocities:
	if any(isnan(vL)) || any(isinf(vL)) || ~any(isreal(vL))
        warning('Node %d of %d has an invalid nodal velocity of: [%e, %e, %e]', ...
            n,numnodes, vL(1),vL(2),vL(3));
	end
end


function [noiseMatrix] = BuildNoiseMatrix(...
    inputMatrix, ...
    mobEps)
	% INPUT:
	%		inputMatrix: square matrix B -- size (sysDim,sysDim)
	%		mobEps: numerical tolerances -- struct
	% OUTPUT:
	%		noiseMatrix: noise matrix -- size (sysDim,sysDim)
    
    % Adding noiseMatrix to B will increase the condition number k(B) and
    % thus the stablility of the algorithm which solves Bv = f, although
    % the accuracy of the solution will be diminished.
    
    % Extraction:
    eps_rCond = mobEps.eps_rcond;
    eps_Noise = mobEps.eps_noise;
    eps_NoiseMax = mobEps.eps_noisemax;
    
    % Check that input matrix has 2 dimensions and is square:
    sysDims = size(inputMatrix); % = [3 3] for full drag matrix
    if size(sysDims,2) ~= 2
        error('Number of system dimensions is %d, but should be %d.', ...
            size(sysDims,2), 2);
    else
        if sysDims(1) ~= sysDims(2)
            error('System size is %d by %d, but should be square.', ...
                sysDims(1), sysDims(2));
        end
    end
    sysDim = sysDims(1); % Size of first dimension of input matrix
    
    targetMatrix = inputMatrix; % Input matrix + noise matrix
    
    % Begin loop to add noise:
    while rcond(targetMatrix) < eps_rCond
        
        % Add uniform isotropic noise:
        targetMatrix = targetMatrix + BuildUniformIsotropicNoiseMatrix(sysDim, mobEps);
        
        if eps_Noise >= eps_NoiseMax
            warning('Poorly conditioned drag matrix; noise added to linear system. Noise parameter: %e', ...
                eps_Noise);
            break
        end
        
        % Increase noise constant until tolerance is met:
        eps_Noise = eps_Noise * 1e1;
    end
    
    % Output only the noise which is added to inputMatrix:
    noiseMatrix = targetMatrix - inputMatrix;
end


function [noiseMatrix] = BuildUniformIsotropicNoiseMatrix(...
    sysDim, ...
    mobEps)
	% INPUT:
    %		sysDim: square linear system dimension -- size (1)
	%		mobEps: numerical tolerances -- struct
	% OUTPUT:
	%		noiseMatrix: matrix of uniform isotropic noise -- size (sysDim,sysDim)
    
    eps_Noise = mobEps.eps_noise;
    noiseMatrix =  eps_Noise * eye(sysDim);
end

