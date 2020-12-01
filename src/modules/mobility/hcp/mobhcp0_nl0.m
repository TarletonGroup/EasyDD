%===============================================================%
% Daniel Hortelano Roig (30/11/2020)
% daniel.hortelanoroig@materials.ox.ac.uk

% Nonlinear HCP mobility law function.
% Solves for v:
% f = (c*|v|^m)*(B*v)

% f: nodal force
% c, m: constants inducing nonlinearity
% B: drag matrix
% v: nodal velocity
% |v|: nodal velocity 2-norm

% For an overview of the HCP mobility law and notation used in this script,
% see Aubry et al: http://dx.doi.org/10.1016/j.jmps.2016.04.019
%===============================================================%

function [vn,fn] = mobhcp0_nl0(fseg,rn,links,connectivity,nodelist,conlist, ...
mobstruct)
% OUTPUT:
%		fn: matrix of nodal forces -- size (N,3) for N input nodes
%       vn: matrix of corresponding nodal response velocities -- size (N,3)

%% Mobility structure extraction %%

%%% Glissile cartesian slip systems:
slipsystems = mobstruct.slipsystemsref; % Cartesian
[numburgs,numplanetypes,~] = size(slipsystems);
maxedgedatadirs = mobstruct.maxedgedatadirs;

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

%%% Nonlinearity:
nlconst = mobstruct.nlconst; % Scaling term
nlexp = mobstruct.nlexp; % Velocity exponent term

%%% MATLAB release:
MATLABrelease = mobstruct.MATLABrelease;

%%% Dimension checks:
assert(all(size(burgscart) == [numburgs,3]), ...
    'Assertion failed: all(size(burgscart) == [numburgs,3])')
assert(all(size(glidecoefficients) == [numburgs,numplanetypes,2]), ...
    'Assertion failed: all(size(glidecoefficients) == [numburgs,numplanetypes,2])')
assert(all(size(crsscoefficients) == [numburgs,numplanetypes + 1]), ...
    'Assertion failed: all(size(crsscoefficients) == [numburgs,numplanetypes + 1])')


%% Mobility law parameters %%

%%% Numerical tolerances:
eps_seglen = 1e7; % Max segment length before linprog could have problems
eps_reflist = 1e-2; % Min difference between actual and reference vectors
eps_decomp = 1e-4; % Min magnitude of decomposed segment component
eps_rconddrag = 1e-9; % Min drag matrix inverse condition number

%%% Nonlinear system analysis:
eps_nrjac = 1e-10; % Min NR Jacobian entries before iteration is reset
eps_nrtol = 1e-14; % Max NR error tolerance allowed in each iteration
nrmaxiter = 20; % Max number of NR iterations before reporting failure
nrlinesearch = 1e-2; % NR line search to improve robustness

%%% Linear system analysis:
eps_dragline = 1/dragline; % Drag line coefficient to improve conditioning
eps_noise = 1e-5; % Noise parameter to improve conditioning
eps_noisemax = 1e-2; % Maximum allowed value for eps_rcondnoise


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
	seglenscrew = zeros(numNbrs,1); % Screw length component
	seglenedges = zeros(numNbrs,maxedgedatadirs); % Edge length components
	unitburg = zeros(numNbrs,3); % Screw data direction
	unitdatadirsedges = zeros(numNbrs,maxedgedatadirs,3); % Edge data directions
	
	% Glide drag coefficients
	screwdragcoeff = zeros(numNbrs,1);
	edgedragcoeff = zeros(numNbrs,maxedgedatadirs);
	
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
        bindex = DetermineIndexOfVectorInReferenceList(burg(i,:),burgscart,eps_reflist);
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
            [glissilenormals,glissileplanetypes] = ...
                CartesianSlipSystemsToPlaneNormalsAndTypes(slipsystems(bindex,:,:));
            
            numglissilenormals = size(glissilenormals,1);
            
            % Match the normal vector to one in the reference list:
            pindex = DetermineIndexOfVectorInReferenceList(ndir(i,:),glissilenormals,eps_reflist);
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
		
        
        %% Determine glide drag coefficients %%
        
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
            
            if abs(seglen(i)) > eps_seglen
                warning('Very large segment length of: %e.\n',seglen(i));
            end
            
            % Decompose segment in the data directions:
			[seglenscrew(i), seglenedgestemp, unitburg(i,:), unitdatadirsedgestemp] = ...
            DecomposeSegmentAsDataDirs(glissilenormals,burgscart(bindex,:),seglenvec(i,:),eps_decomp,MATLABrelease);
			
			% Account for padding due to number of data directions available for a segment:
			seglenedges(i,1:size(glissilenormals,1)) = seglenedgestemp;
			unitdatadirsedges(i,1:size(glissilenormals,1),:) = reshape(unitdatadirsedgestemp,1,size(glissilenormals,1),3);
			
			% Screw drag coefficients:
			screwdragcoeff(i) = glidecoefficients(bindex,planetype,1);
			
			% Edge drag coefficients:
            for p = 1:size(glissilenormals,1)
				planetypeedge = glissileplanetypes(p);
				edgedragcoeff(i,p) = glidecoefficients(bindex,planetypeedge,2);
            end
        end
        
        
        %% Determine nodal forces %%
        
        % Contribution to total force on node n0 from segment i:
		fsegn0 = fseg(linkid, 3 * (posinlink - 1) + (1:3));
                
        % Determine CRSS:
        if planetype == 6
            crss = crsssessile; % Sessile CRSS
        else
            crss = crsscoefficients(bindex,planetype);
        end
        % Impose CRSS:
        fsegn0plane = fsegn0 - dot(fsegn0,ndir(i,:)) * ndir(i,:);
        if norm(fsegn0plane) / seglen(i) < crss
            % Segment force was not sufficient to overcome CRSS.
            crsspassed(i) = 0;
            % Treat this segment as a junction:
            junctionexists(i) = 1;
        else
            % Segment force was sufficient to overcome CRSS.
            crsspassed(i) = 1;
        end
        
        % Superpose fsegn0 to the total force on n0:
        fn(n,:) = fn(n,:) + fsegn0;
        
                
	end % Repeat for all connected segments
    
    
    %% Construct total drag matrix %%
    
    Btotal = zeros(3,3);
    
    for i = 1:numNbrs
        
        % Segments with zero length do not contribute:
        if seglen(i) == 0
            continue
        end
        
        % Junction contribution:
        if junctionexists(i) == 1
            
            % Construct the junction drag matrix:
            junctdragmatrix = BuildJunctionDragMatrix(linedir(i,:), dragline, dragsessile);
            
            % Add contribution to total drag function:
            Btotal = Btotal + CalculateDragMatrixContribution(seglen(i),junctdragmatrix);
            
            continue
        end

        % Screw contribution:
        if abs(seglenscrew(i)) > 0
            
            % Unit direction vectors:
            linevector = unitburg(i,:);
            climbvector = ndir(i,:);
            glidevector = cross(linevector,ndir(i,:))/norm(cross(linevector,ndir(i,:)));
            
            % Construct the screw drag matrix:
            screwdragmatrix = BuildGlissileDragMatrix(linevector, climbvector, glidevector, ...
                                                   dragline, dragclimb, screwdragcoeff(i));
            
            Btotal = Btotal + CalculateDragMatrixContribution(seglenscrew(i),screwdragmatrix);
        end
        
        % Edge contribution:
        for p = find(seglenedges(i,:)) % Cycles through indices corresponding to nonzero entries
            
            % Unit direction vectors:
            linevector = reshape(unitdatadirsedges(i,p,:),1,3);
            climbvector = cross(unitburg(i,:),linevector)/norm(cross(unitburg(i,:),linevector));
            glidevector = unitburg(i,:);
            
            % Construct the edge drag matrix:
            edgedragmatrix = BuildGlissileDragMatrix(linevector, climbvector, glidevector, ...
                                                  dragline, dragclimb, edgedragcoeff(i,p));
            
            Btotal = Btotal + CalculateDragMatrixContribution(seglenedges(i,p),edgedragmatrix);
        end
    end
    
    
    %% Solve nonlinear drag system
    % Performs Newton-Raphson iteration.
    
    % Iteration: J(k)*(v(k+1) - v(k)) = -F(k)
	% rootfn = F(k)             	<-- Root function
	% jacfn = J(k) = grad(F(k))     <-- Jacobian of root function
    % nrerror2 = |rootfn|^2         <-- Squared residual norm (error)
    
    % Initial (linear) guess of solution:
    vNR = (Btotal \ fn(n,:)')'; % Velocity (target)
    rootfn = RootFunction(fn(n,:), nlconst, nlexp, vNR, Btotal); % Root function
    jacfn = JacobianFunction(nlconst, nlexp, vNR, Btotal); % Jacobian function
    nrerror2 = norm(rootfn)^2;
    
    nrcount = 0;
    % Do NR iteration until tolerance met or max iterations reached:
    while (nrcount < nrmaxiter) && (nrerror2 > eps_nrtol)
        
        % Try to resolve potential singularity induced by stationary point:
        if  all(abs(jacfn) < eps_nrjac,'all')            
            vNR = 10*(rand(1,3)-0.5); % Try random velocity
            rootfn = RootFunction(fn(n,:), nlconst, nlexp, vNR, Btotal);
            jacfn = JacobianFunction(nlconst, nlexp, vNR, Btotal);
        end
        
        % Solve nonlinear system:
        deltav = (jacfn \ (-rootfn'))';
        vNR = vNR + nrlinesearch * deltav;
        
        % Compute corresponding root and jacobian functions:
        rootfn = RootFunction(fn(n,:), nlconst, nlexp, vNR, Btotal);
        jacfn = JacobianFunction(nlconst, nlexp, vNR, Btotal);
        
        nrerror2 = norm(rootfn)^2;
        
        nrcount = nrcount + 1;
    end
    
    
    %% Analyse nonlinear drag system
    
    if (nrcount >= nrmaxiter) || (nrerror2 > eps_nrtol)
        % Newton-Raphson iteration cannot converge for node n.
        % The system will instead be approximated as linear.
        
        warning('Newton-Raphson has failed for node: %d of %d. The system will be approximated as linear.\n', ...
            n, numnodes);
        fprintf('nrcount: %d || det(Btotal): %0.16s || nrerror2: %0.16s || nrerror2tol: %0.16s\n', ...
            nrcount, det(Btotal), nrerror2, eps_nrtol);
    else
        % Newton-Raphson iteration converged for node n.
        
        %%%%%%%%%%%%%%
        vn(n,:) = vNR; % Nonlinear system
        %%%%%%%%%%%%%%
        
        continue % Continue to next node
    end
    
    
    %% Analyse and solve linear drag system %%
    
	%%% Verify properties of drag matrix:
    if ~isreal(Btotal)
        warning('The drag matrix should be real, but is not.');
    end
    if ~issymmetric(Btotal)
        warning('The drag matrix should be symmetric, but is not.');
    end
    
    %%% Prepare for diagnosis and potential restructuring of linear system:
    % Diagnostic information:
    [evecs,evals] = eigs(Btotal);
    imageidxs = (1:size(evals,1))'; % Default (= [1 2 3]')
    % Restructured matrices:
    BtotalF = Btotal;
    fnF = fn(n,:);
    
    %%% Analyse linear system singularities:
    if any(diag(evals) < eps(10))
        % Drag matrix contains at least one zero eigenvalue to machine
        % precision; corresponding eigenvectors are nullspace basis
        % vectors. The linear system will be solved in its eigenbasis
        % with reduced dimensionality (nullspace dimensions removed) to
        % improve stability.
        
        % NOTE: this disables nodal motion in nullspace direction vectors.
        
        % Set near-zero eigenvalues to zero:
        evals(abs(evals) < eps(10)) = 0; % eps(10) \approx 1.78e-15
        
        imageidxs = diag(evals) ~= 0; % Nonzero eigenvalue indices
        
        % Rotate nodal force vector into eigenbasis:
        fnF = (evecs * fnF')';
        
        % Remove dimensions corresponding to zero eigenvalues:
        fnF = fnF(imageidxs);
        BtotalF = evals(imageidxs,imageidxs);
    end
    dimsys = size(Btotal,1); % Original system dimension
    dimsysF = size(BtotalF,1); % Restructured system dimension
    
    %%% Analyse drag matrix conditioning:
    if rcond(BtotalF) < eps_rconddrag
        % Drag matrix is poorly conditioned; try decreasing the disparity
        % of the magnitude of its set of eigenvalues.
        
        % Increase drag in the line directions of the connected segments:
        for i = 1:numNbrs
            
            seglinevec = linedir(i,:);
            if dimsysF < dimsys
                
                % Rotate line vector into eigenbasis:
                seglinevec = (evecs * seglinevec')';
                % Reduce dimension:
                seglinevec = seglinevec(imageidxs);
            end
            BtotalF = BtotalF + BuildJunctionDragMatrix(seglinevec,eps_dragline,0);
        end
        
        if rcond(BtotalF) < eps_rconddrag
            % Linear system is still poorly conditioned; increase the
            % stiffness of the node in isotropic directions by adding
            % unform isotropic noise (i.e. scaled identity matrix).
            
            BtotalF = BtotalF + BuldUniformIsotropicNoiseMatrix(BtotalF, ...
                eps_rconddrag,eps_noise,eps_noisemax);
        end
    end
    
    
    %%% Final nodal velocity
    
    vnF = (BtotalF\fnF')';
    
    if dimsysF < dimsys
        
        % Replace the nullspace dimension with zeros:
        vnF0 = zeros(1,3);
        vnF0(imageidxs) = vnF;
        % Undo the rotation into the eigenbasis:
        vnF0 = (evecs' * vnF0')'; % Note orthogonality of evecs
        
        vnF = vnF0;
    end
    
    % Check consistency of system size:
    assert(size(Btotal,1) == size(vnF',1), ...
        'Assertion failed: size(Btotal,1) == size(vnF,2)')
    assert(size(BtotalF,1) == size(fnF',1), ...
            'Assertion failed: size(BtotalF,1) == size(fnF,2)')
    assert(size(imageidxs,1) == size(fnF',1), ...
            'Assertion failed: size(imageidxs,1) == size(fnF,2)')
    
    % Check validity of nodal velocities:
	if any(isnan(vnF)) || any(isinf(vnF)) || ~any(isreal(vnF))
        warning('Node %d of %d has an invalid nodal velocity of: [%e, %e, %e]', ...
            n,numnodes, vnF(1),vnF(2),vnF(3));
	end
    
    %%%%%%%%%%%%%%
    vn(n,:) = vnF; % Linear system
    %%%%%%%%%%%%%%
    
    
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

function [vIndex] = DetermineIndexOfVectorInReferenceList(cartVec,cartRefVecs,eps_RefList)
	% INPUT:
	%		cartVec: input vector -- size (1,3)
	% 		cartRefVecs: matrix of vectors in reference list -- size (P,3)
    %       eps_RefList: scalar tolerance for whether cartVec is in reference list -- size (1)
	% OUTPUT:
	%		vIndex: vector of indexes of all matched vectors in reference list -- size (1,Q)
    
    % NOTE: if no matches are found, vIndex = [].
    
    epsVec = eps_RefList;
    
    numRefVecs = size(cartRefVecs,1);
    unitCartVec = cartVec/norm(cartVec);
    
    vIndex = zeros(1,numRefVecs); % Preallocate for max possible matches
    for v = 1:numRefVecs
        unitCartRefVec = cartRefVecs(v,:)/norm(cartRefVecs(v,:));
        if norm(cross(unitCartVec,unitCartRefVec)) < epsVec
            
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
	% MB indices vertically concatenated in the same order as those hard-coded in slipsSystems:
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

function [segLenDecompScrew, segLenDecompEdge, unitBurgVec, unitDataDirsMatrixEdges] = ...
    DecomposeSegmentAsDataDirs(glissilenormals,burgVec,segLenVec,eps_Decomp,MATLABRelease)
	% INPUT:
	%		glissilenormals: glissile plane normal vectors for the current segment -- size (P,3) for P glissile planes
	% 		burgVec: Burgers vector of the segment -- size (1,3)
	% 		segLenVec: length vector of the segment -- size (1,3)
    %       eps_Decomp: numerical tolerance for length component of a decomposed segment  -- size (1)
    %       MATLABRelease: string of current running MATLAB client release version (e.g. '(R2020a)')
	% OUTPUT:
	%		segLenDecompScrew: screw component of the segment length vector decomposed in the data directions -- size (1)
	%		segLenDecompEdge: edge components "													   " -- size (1,P)
	%		unitBurgVec: screw data direction (also the unit Burgers vector of the current segment) -- size (1,3)
	%		unitDataDirsMatrixEdges: edge data directions -- size (P,3)
		
	numEdgesDataDirs = size(glissilenormals,1);
	numDataDirs = numEdgesDataDirs + 1;
	unitDataDirsMatrix = zeros(3,numDataDirs);
	
	% Screw direction:
	unitDataDirsMatrix(:,1) = (burgVec/norm(burgVec))';
	
	% Edge directions:
	for k = 1:numEdgesDataDirs
		normal = glissilenormals(k,:);
		dirVec = cross(burgVec,normal);
		dirUnitVec = dirVec/norm(dirVec);
		unitDataDirsMatrix(:,k+1) = dirUnitVec';
	end
	segLenDecomp = LinearProgrammingL1Minimisation(unitDataDirsMatrix,segLenVec',MATLABRelease);
	
	% Set negligible segment decomposition components to zero:
	segLenDecomp(abs(segLenDecomp)/norm(segLenVec) < eps_Decomp) = 0;
	
	% Separate screw and edge components for output:
	segLenDecompScrew = segLenDecomp(1)';
	segLenDecompEdge = segLenDecomp(2:end)';
	unitBurgVec = unitDataDirsMatrix(:,1)';
	unitDataDirsMatrixEdges = unitDataDirsMatrix(:,2:end)';
end


 function [outVector] = ...
     LinearProgrammingL1Minimisation(matrixOfDirections,inVector,MATLABRelease)
	% INPUT:
	%		matrixOfDirections: matrix of direction vectors -- size (dimS,numDD) for numDD direction vectors
	% 		inVector: vector to decompose -- size (dimS,1)
    %       MATLABRelease: string of current running MATLAB client release version (e.g. '(R2020a)')
	% OUTPUT: 
	%		outVector: approximately the same vector as inVector, but decomposed as a subset of the
    %                  data direction vectors in matrixOfDirections -- size (numDD,1)
	
	% Optimization condition: minimise the L1 norm of outVector.
    
    % After some derivation, the decomposition problem may be formulated as
    % a linear program with linear constraints.
    % When properly formulated, one can solve this using MATLABs linprog:
    % https://uk.mathworks.com/help/optim/ug/linprog.html
    
    
    %%% Build linprog arguments (for L1 norm optimization)
    
    % Number of spatial dimensions:
    dimS = size(matrixOfDirections,1);
    % Number of data directions:
    numDD = size(matrixOfDirections,2);
    
    % Initialize vectors for linprog:
	vectorOfOnes = ones(numDD,1);
    vectorOfZeros = zeros(numDD,1);
    matrixOfZeros = zeros(dimS,numDD);
    
    % Constraint condition:
    constraint = [vectorOfZeros; vectorOfOnes]';
	
    % Equality condition:
	Aeq = [matrixOfDirections matrixOfZeros];
    beq = inVector;
    
    % Inequality condition:
    Aneq = [eye(numDD) -eye(numDD); -eye(numDD) -eye(numDD)];
    bneq = [vectorOfZeros; vectorOfZeros]';
    
    % outVector lower and upper limits:
    lb = [-Inf * vectorOfOnes; vectorOfZeros]; % Lower limit due to L1 norm minimisation constraint
    ub = [Inf * vectorOfOnes; Inf * vectorOfOnes]; % No upper limit
    
    
    %%% Set up linprog options
    
    % NOTE: linprog has a bug (error code -2_4) when using its default
    % dual-simplex algorithm, which was fixed as of MATLAB R2019a.
    % One workaround (officially dev-approved) is to use the interior-point
    % algorithm instead of dual-simplex.
    MATLABReleaseYear = str2double(MATLABRelease(3:6)); % e.g. '2020'
    
    if MATLABReleaseYear < 2019
        
        % Use interior-point algorithm for bug workaround:
        options = optimoptions('linprog','Algorithm','interior-point');
        
        % Tolerances:
        eps_OptTol = 1e-6; % Default OptimalityTolerance
        eps_OptTolMax = 1e-1; % Maximum allowed OptimalityTolerance
        
    elseif MATLABReleaseYear >= 2019
        
        % Use the dual-simplex algorithm:
        options = optimoptions('linprog','Algorithm','dual-simplex');
        
        % Tolerances:
        eps_OptTol = 1e-7; % Default
        eps_OptTolMax = 1e-2;
    end
    
    % Generic options for both linprog algorithms:
    eps_MaxOptdT = 1.0e0; % Tolerance for max linprog computation time
    options.Display = 'off'; % Suppresses optimization-related output
    options.OptimalityTolerance = eps_OptTol;
    %options.MaxIterations = 200; % Default max linprog iterations
    %options.MaxTime = Inf; % Default max linprog computation time
    
    % Solve optimisation problem:
    preOptT = tic;
	[optimisedVec,fVal,exitFlag,outputLin,lambdaLin] = linprog(constraint,Aneq,bneq,Aeq,beq,lb,ub,options);
    postOptT = toc;
    optdT = postOptT - preOptT;
    % NOTE: optimizedVec has size (2 * numDD,1). The entries are:
    % optimizedVec(numDD:2*numDD) = abs(optimizedVec(1:numDD))
    
    if optdT > eps_MaxOptdT
        warning('Large linprog computation time of: %e seconds', ...
            optdT);
    end
    
    
    %%% Check linprog convergence
        
    while exitFlag <= 0
        % linprog failed to converge within tolerances.
        % Run linprog with less restrictive tolerances until successful
        % convergence or failure reached beyond maximum allowed tolerance.
        
        if eps_OptTol >= eps_OptTolMax
            
            % linprog diagnostic information:
            fprintf('eps_OptTolMax: '); disp(eps_OptTolMax);
            fprintf('\neps_OptTol: '); disp(eps_OptTol);
            fprintf('\nlambda: '); disp(lambdaLin);
            fprintf('\nfVal: '); disp(fVal);
            fprintf('\n=== output ===\n'); disp(outputLin);
            
            error('linprog attempted to decompose the vector: [%e, %e, %e], but failed within tolerances. Exit flag: %d.', ...
            inVector(1),inVector(2),inVector(3), exitFlag);
        end

        eps_OptTol = eps_OptTol * 1e1;
        options.OptimalityTolerance = eps_OptTol;
        
        % linprog with new tolerances while tracking computation time:
        preOptTConv = tic;
        [optimisedVec,fVal,exitFlag,outputLin,lambdaLin] = linprog(constraint,Aneq,bneq,Aeq,beq,lb,ub,options);
        postOptTConv = toc;
        optdTConv = postOptTConv - preOptTConv;
        
        % If linprog computation time is too large, reduce MaxTime:
        if optdTConv > eps_MaxOptdT
            options.MaxTime = eps_MaxOptdT;
        end
    end
    
    % Successful convergence!
    outVector = optimisedVec(1:numDD);
end


%% Build drag matrix %%

function [dragMatrixContribution] = CalculateDragMatrixContribution(segLength,dragMatrix)
	
	dragMatrixContribution = 0.5 * abs(segLength) .* dragMatrix;
end


function [junctDragMatrix] = BuildJunctionDragMatrix(lineVector, junctLineDragCoeff, junctGlideDragCoeff)
	
    sysDim = size(lineVector,2); % = 3 for full linear system
    
	junctLineDirMatrix  = lineVector'  * lineVector;
	junctGlideDirMatrix = eye(sysDim) - junctLineDirMatrix;
	
	junctDragMatrix = junctLineDragCoeff * junctLineDirMatrix + junctGlideDragCoeff * junctGlideDirMatrix;
end


function [dragMatrix] = BuildGlissileDragMatrix(lineVector, climbVector, glideVector, lineDragCoefficient, climbDragCoefficient, glideDragCoefficient)

	lineDirMatrix = lineVector'  * lineVector;
	climbDirMatrix = climbVector' * climbVector;
	glideDirMatrix = glideVector' * glideVector;
	
	dragMatrix = lineDragCoefficient * lineDirMatrix + climbDragCoefficient * climbDirMatrix + glideDragCoefficient * glideDirMatrix;
end


%% Build Jacobian matrix %%

function [rootFn] = RootFunction(f,c,m,vel,dragMatrix)
    
    % Root function:
    % F = f - (c*|v|^m)*(B*v)
    
    % f: size (1,3)
    % c,m: size (1)
    % vel: size (1,3)
    % dragMatrix: size (3,3)
    
    %%%
    
    % Change to size (3,1):
    f = f';
    vel = vel';
    
    normvel = norm(vel);
    
    rootFn = f - c*(normvel^(m))*(dragMatrix*vel);
    
    rootFn = rootFn'; % Change to size (1,3)
end


function [jacFn] = JacobianFunction(c,m,vel,dragMatrix)
    
    % Jacobian function:
    % J = -c*|v|^(m)*B - c*m*|v|^(m-2)*(B*v)*v
    
    % c,m: size (1)
    % vel: size (1,3)
    % dragMatrix: size (3,3)
    
    % Warning: stationary point at |v| = 0 induces singularity in NR
    
    %%%
    
    % Change to size (3,1):
    vel = vel';
    
    normvel = norm(vel);
    
    jacFn = -c*(normvel^(m))*dragMatrix - c*m*(normvel^(m-2))*(dragMatrix*vel)*vel';
    
    jacFn = jacFn'; % Change to size (1,3)
end


%% Analyse linear system %%

function [noiseMatrix] = BuldUniformIsotropicNoiseMatrix(inputMatrix,eps_rCond,eps_Noise,eps_NoiseMax)
	% INPUT:
	%		inputMatrix: square matrix B -- size (sysDim,sysDim)
	%		eps_rCond: minimum value for rcond(inputMatrix) -- size (1)
    %		eps_Noise: noise parameter -- size (1)
    %		eps_NoiseMax: maximum allowed noise parameter -- size (1)
	% OUTPUT:
	%		noiseMatrix: matrix of uniform isotropic noise -- size (sysDim,sysDim)
    
    % Adding noiseMatrix to B will increase the condition number k(B) and
    % thus the stablility of the algorithm which solves Bv = f, although
    % the accuracy of the solution will be diminished.
    
    sysDim = size(inputMatrix,1); % = 3 for full drag matrix
    targetMatrix = inputMatrix;
    
    while rcond(targetMatrix) < eps_rCond
        
        % Add noise to input matrix:
        targetMatrix = targetMatrix + eps_Noise * eye(sysDim);
        
        if eps_Noise >= eps_NoiseMax
            warning('Poorly conditioned drag matrix; isotropic stiffness added to linear system. Noise parameter: %e', ...
                eps_Noise);
            break
        end
        
        % Increase noise constant until tolerance is met:
        eps_Noise = eps_Noise * 1e1;
    end
    
    % Output only the noise which is added to inputMatrix:
    noiseMatrix = targetMatrix - inputMatrix;
end

