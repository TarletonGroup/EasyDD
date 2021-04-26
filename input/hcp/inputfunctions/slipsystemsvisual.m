%===============================================================%
% Daniel Hortelano Roig (01/02/2020)
% daniel.hortelanoroig@materials.ox.ac.uk

% This script is used for visualising the slip systems of HCP metals.
% Visualises zirconium (Zr) by default.

% To select the specific slip system to visualise,
% select bidx (Burgers vector index) and/or pidx (normal vector index).
%===============================================================%

%% Determine slip systems

% Define slip system parameters:
numburgs = 10;
numplanetypes = 5;
HCPa = 1.0;
HCPc = 1.593; % Zr
refHCPa1 = HCPa * [-1/2 sqrt(3)/2 0];
refHCPa2 = HCPa * [-1/2 -sqrt(3)/2 0];
refHCPa3 = HCPa * [1 0 0];
refHCPa4 = HCPc * [0 0 1];
refHCPmat = [refHCPa1' refHCPa2' refHCPa4']; % Linearly independent columns

% Define slip systems:
slipsystems = SlipSystemsDef(numburgs, numplanetypes);
slipsystemscart = SlipSystemsToCartesian(slipsystems, refHCPa1,refHCPa2,refHCPa3,refHCPa4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select desired slip system to visualise:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bidx = 7; % Burgers vector index
pidx = 4; % Plane index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select crystal axes rotation (rotates all vectors):
ssangx = 0/180 * pi; % Rotation angle
ssRx = [1 0 0; 0 cos(ssangx) -sin(ssangx); 0 sin(ssangx) cos(ssangx)]; % Rotation matrix (wrt x as example)
rotcrystalmat = ssRx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Compute associated vectors and points

burgmb = slipsystems{bidx,pidx,1};
planemb = slipsystems{bidx,pidx,2}(1,:);
planemill = planemb; planemill(3) = [];

burgcart = slipsystemscart{bidx,pidx,1};
planecart = slipsystemscart{bidx,pidx,2}(1,:);

planepts = [refHCPmat(:,1)/planemill(1), refHCPmat(:,2)/planemill(2), refHCPmat(:,3)/planemill(3)];

% Set non-axes-intersecting plane points to very far away:
eps_dist = 1e10;
infidx = [0 0 0]; % =1: set to pseudoinfinity, =0: don't
for i = 1:size(planepts,2)
    if any(isinf(planepts(:,i))) || any(isnan(planepts(:,i)))
        planepts(:,i) = refHCPmat(:,i)*eps_dist;
        infidx(i) = 1;
    end
end

% Set plane "mid point":
planemidpt = [0 0 0]'; % Initialisation
for i = 1:size(planepts,2)
    if infidx(i) == 0 % Only superpose non-pseudoinfinite vertices
        planemidpt = planemidpt + planepts(:,i); % Attain centroid
    end
end
numnoninf = sum(infidx == 0); % Number of non-pseudoinfinitely-far-away vertices
planemidpt = planemidpt / numnoninf;

% Rotating all vectors:
planepts = rotcrystalmat * planepts;
planemidpt = rotcrystalmat * planemidpt;
planecart = (rotcrystalmat * planecart')';
burgcart = (rotcrystalmat * burgcart')';
refHCPa1 = (rotcrystalmat * refHCPa1')';
refHCPa2 = (rotcrystalmat * refHCPa2')';
refHCPa3 = (rotcrystalmat * refHCPa3')';
refHCPa4 = (rotcrystalmat * refHCPa4')';


%% Plot figure

%close all
% Setup figure:
%figname = append('Slip system: (b',num2str(bidx),')[ ',num2str(burgmb,'%1.2f'),' ] (p',num2str(pidx),')[ ',num2str(planemb,'%0.0f'),' ]');
%figname = append('Slip system: (b',num2str(bidx),')[ ',compose("%0.2f",burgmb),' ] (p',num2str(pidx),')[ ',compose("%0.2f",planemb),' ]');
figname = append('Slip system: b',num2str(bidx),'[ ',rats(burgmb,5),' ]    p',num2str(pidx),'( ',rats(planemb,5),' )');
figure('Name',figname);
grid on, hold on

% xyz axes:
b1 = plot3([-3 3], zeros(2),zeros(2), ... % x-axis
    '--','LineWidth',0.3,'Color','black');
plot3(zeros(2), [-3 3], zeros(2), ... % y-axis
    '--','LineWidth',0.3,'Color','black');
plot3(zeros(2), zeros(2), [-3 3], ... % z-axis
    '--','LineWidth',0.3,'Color','black');

p1 = fill3(planepts(1,:),planepts(2,:),planepts(3,:),'r'); % Plane
alpha(0.3)
q1 = quiver3(planemidpt(1),planemidpt(2),planemidpt(3),burgcart(1),burgcart(2),burgcart(3), ...
        'Color','blue','LineWidth',1, ...
        'ShowArrowHead','on','MaxHeadSize',1e1); % Burgers vector
q2 = quiver3(planemidpt(1),planemidpt(2),planemidpt(3),planecart(1),planecart(2),planecart(3), ...
        'Color','red','LineWidth',1, ...
        'ShowArrowHead','on','MaxHeadSize',1e1); % Plane normal vector

% HCP basis vectors:
b2 = quiver3(0,0,0,refHCPa1(1),refHCPa1(2),refHCPa1(3), ...
        '--','Color','magenta','LineWidth',2, ...
        'ShowArrowHead','off','MaxHeadSize',1e0);
quiver3(0,0,0,refHCPa2(1),refHCPa2(2),refHCPa2(3), ...
        '--','Color','magenta','LineWidth',2, ...
        'ShowArrowHead','off','MaxHeadSize',1e0);
quiver3(0,0,0,refHCPa3(1),refHCPa3(2),refHCPa3(3), ...
        '--','Color','magenta','LineWidth',2, ...
        'ShowArrowHead','off','MaxHeadSize',1e0);
quiver3(0,0,0,refHCPa4(1),refHCPa4(2),refHCPa4(3), ...
        '--','Color','magenta','LineWidth',2, ...
        'ShowArrowHead','off','MaxHeadSize',1e0);

%%% HCP hexagonal prism

% Generating vertices:
nH = 6;
zl = 0; hH = HCPc;
AH = ones(nH + 1);
z1H = AH(:,1) * zl;
z2H = z1H + hH;
tH = 0 : 2*pi/nH : 2*pi;
xH = HCPa * cos(tH);
yH = HCPa * sin(tH);

% Rotating vertices:
HCPbot = rotcrystalmat * [xH; yH; z1H']; % Rotate vertices on bottom of hexagonal prism
HCPtop = rotcrystalmat * [xH; yH; z2H']; % Rotate vertices on top of hexagonal prism
xHbot = HCPbot(1,:); xHtop = HCPtop(1,:);
yHbot = HCPbot(2,:); yHtop = HCPtop(2,:);
z1Hbot = HCPbot(3,:)'; z2Htop = HCPtop(3,:)'; % New vertex z-coordinates

% Plotting rotated hexagonal prism:
surf([xHbot;xHtop].', [yHbot;yHtop].', [z1Hbot,z2Htop], ...
    'FaceColor','cyan','EdgeColor',[0.4940 0.1840 0.5560],'LineStyle','-.','FaceAlpha',0.2) % Fill hexagonal prism sides
patch('XData',[xHbot';xHtop'], 'YData',[yHbot';yHtop'], 'ZData',[z1Hbot,z2Htop], ...
    'FaceColor','cyan','EdgeColor',[0.4940 0.1840 0.5560],'LineStyle','-.','FaceAlpha',0.2) % Fill polygon bot and top

%%% Axes and labels

xlim([-3 3]);
ylim([-3 3]);
zlim([-3 3]);
xlabel('x','FontSize',14); % Label axes
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
title(figname,'FontSize',14); % Plot title

lgd = legend([b2 p1 q2 q1],{'HCP axes', 'Plane', 'Plane normal', 'Burgers vector'}, ...
    'Location','southeast'); % Legend
lgd.FontSize = 12;
lgd.Title.String = 'Legend of Vectors';

axis('square')



%% Function definitions

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

function [planeNormals] = PlanesMBToCartesianNormals(planesRefmb,a1,a2,~,a4)
	% INPUT:
	%       planesRefmb: matrix of MB indices (hkil) describing planes -- size (Q,4) for Q slip systems
	%                     NOTE: a size (Q,3) matrix of (hkl) Miller indices will also work
	%       a1,a2,a4: HCP lattice vectors -- each size (1,3)
    %                     NOTE: a1,a2,a4 (a3 omitted) must be linearly dependent!
	% OUTPUT: 
	%		planeNormals: matrix of unit normal vectors of the planes generated by the corresponding MB indices -- size (Q,3)
	
	planeNormals = zeros(size(planesRefmb,1),3);
	
	for i = 1:size(planesRefmb,1)
		
		indices = planesRefmb(i,:); % MB indices describing a plane in HCP
		
		% Notice that (a1,a2,a4) is the chosen linearly independent HCP basis:
		unitNormalVec = IndicesMBToCartesianNormal(indices,a1,a2,a4);
		
		planeNormals(i,:) = unitNormalVec; % Unit normal vector in cartesian coordinates
	end
end

function [unitNormalVec] = IndicesMBToCartesianNormal(indices,a1,a2,a3)
	% INPUT:
	%		indices: four MB indices (hkil) describing the plane -- size (1,4)
	%                NOTE: Miller indices (hkl) with size (1,3) will also work
	%       a1,a2,a3: three linearly independent HCP lattice basis vectors -- each size (1,3)
    %       eps_VolLat: tolerance for linear independence of the three input HCP lattice vectors -- size (1)
	% OUTPUT:
	%		unitNormalVec: corresponding plane unit normal vector in cartesian coordinates -- size (1,3)
    
    % Tolerance for linear independence of three input HCP vectors:
    eps_VolLat = 1e-12;
    
	% HCP lattice basis vectors in a matrix form:
	basis = [a1;a2;a3]; % size (3,3)
	
	% Volume of cell generated by (a1,a2,a3):
	vol = det(basis); % size (1)
	
	% Checking linear independence of (a1,a2,a3):
	if abs(vol) < eps_VolLat
        warning('The volume generated by the three chosen HCP lattice vectors is %0.2s, which is less than the tolerance %0.2s.\n',abs(vol),eps_VolLat)
	end
	
	% Reciprocal lattice basis vectors:
	g1 = cross(a2,a3) ./ vol;
	g2 = cross(a3,a1) ./ vol;
	g3 = cross(a1,a2) ./ vol;
	
	if size(indices,2) == 4 % indices contains Miller-Bravais indices
		if indices(3) ~= -(indices(1) + indices(2))
            warning('i = %d != -(h + k) = %d for these Miller-Bravais indices.\n',indices(3),-(indices(1) + indices(2)))
		end
		indices(3) = []; % Removes linearly dependent index
	elseif size(indices,2) == 3 % indices contains Miller indices
		% Do nothing! Miller indices are apparently already given
    else
        error('The number of indices describing this plane is unexpectedly %d, and should be either 3 or 4.\n',size(indices,2))
	end
	
	% Miller indices in standard notation:
	h = indices(1);
	k = indices(2);
	l = indices(3);
	
	% Normal vector to (hkl) plane:
	normalVec = h*g1 + k*g2 + l*g3;
	
	% Unit normal vector:
	unitNormalVec = normalVec/norm(normalVec);
end