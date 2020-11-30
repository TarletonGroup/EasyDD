%===============================================================%
% Daniel Celis Garza & Daniel Hortelano Roig (29/11/2020)

% Ensure we have all input variables. Also organises essential
% input data into structures.
%===============================================================%

%% Verify input variables

if ~exist('MU', 'var')
    MU = 1;
    fprintf('Providing default value for MU = %f.\n', MU)
end

if ~exist('NU', 'var')
    NU = 0.305;
    fprintf('Providing default value for NU = %f.\n', NU)
end

if ~exist('amag', 'var')
    amag = 3.18e-4;
    fprintf('Providing default value for amag = %f.\n', amag)
end

if ~exist('mumag', 'var')
    mumag = 1.6e5; % MPa only for plotting.
    fprintf('Providing default value for mumag = %f.\n', mumag)
end

if ~exist('dx', 'var')
    dx = 30 / amag; % 30 micron
    fprintf('Providing default value for dx = %f.\n', dx)
end

if ~exist('dy', 'var')
    dy = 5 / amag; % 5 micron
    fprintf('Providing default value for dy = %f.\n', dy)
end

if ~exist('dz', 'var')
    dz = 5 / amag; % 5 micron
    fprintf('Providing default value for dz = %f.\n', dz)
end

if ~exist('mx', 'var')
    mx = 20;
    fprintf('Providing default value for mx = %f.\n', mx)
end

if ~exist('Bcoeff', 'var')
    mobstruct = struct('screw', 10, 'edge', 1, 'climb', 1e10, 'line', 1e-4);
    fprintf('Providing default value for Bcoeff.screw = %f.\n', mobstruct.screw)
    fprintf('Providing default value for Bcoeff.edge = %f.\n', mobstruct.edge)
    fprintf('Providing default value for Bcoeff.climb = %f.\n', mobstruct.climb)
    fprintf('Providing default value for Bcoeff.line = %f.\n', mobstruct.line)
end

if ~exist('maxconnections', 'var')
    maxconnections = 4;
    fprintf('Providing default value for maxconnections = %d.\n', maxconnections)
end

if ~exist('lmax', 'var')
    lmax = 0.25 / amag;
    fprintf('Providing default value for lmax = %f.\n', lmax)
end

if ~exist('lmin', 'var')
    lmin = 0.1 / amag;
    fprintf('Providing default value for lmin = %f.\n', lmin)
end

if ~exist('areamin', 'var')
    areamin = lmin * lmin * sin(60/180 * pi) * 0.5;
    fprintf('Providing default value for areamin = %f.\n', areamin)
end

if ~exist('areamax', 'var')
    areamax = 20 * areamin;
    fprintf('Providing default value for areamax = %f.\n', areamax)
end

if ~exist('CUDA_flag', 'var')
    CUDA_flag = false;
    fprintf('Providing default value for CUDA_flag = %d.\n', CUDA_flag)
end

if ~exist('a_trac', 'var')
    a_trac = true;
    fprintf('Providing default value for a_trac = %d.\n', a_trac)
end

if ~exist('calculateR_hat', 'var')
    calculateR_hat = 1;
    fprintf('Providing default value for calculateR_hat = %d.\n', calculateR_hat)
end

if ~exist('doremesh', 'var')
    doremesh = 1; %flat set to 0 or 1 that turns the remesh functions off or on
    fprintf('Providing default value for doremesh = %d.\n', doremesh)
end

if ~exist('docollision', 'var')
    docollision = 1; %flat set to 0 or 1 that turns collision detection off or on
    fprintf('Providing default value for docollision = %d.\n', docollision)
end

if ~exist('doseparation', 'var')
    doseparation = 1; %flat set to 0 or 1 that turns splitting algorithm for highly connected node off or on
    fprintf('Providing default value for doseparation = %d.\n', doseparation)
end

if ~exist('dovirtmesh', 'var')
    dovirtmesh = 1; %flat set to 0 or 1 that turns remeshing of virtual nodes off or on
    fprintf('Providing default value for dovirtmesh = %d.\n', dovirtmesh)
end

if ~exist('dt0', 'var')
    dt0 = 1e6;
    fprintf('Providing default value for dt0 = %f.\n', dt0)
end

if ~exist('dt', 'var')
    dt = dt0;
    fprintf('Providing default value for dt = %f.\n', dt)
end

if ~exist('simTime', 'var')
    simTime = 0;
    fprintf('Providing default value for simTime = %f.\n', simTime)
end

if ~exist('totalSimTime', 'var')
    totalSimTime = 1e9;
    fprintf('Providing default value for totalSimTime = %f.\n', totalSimTime)
end

if ~exist('curstep', 'var')
    curstep = 0;
    fprintf('Providing default value for curstep = %d.\n', curstep)
end

if ~exist('a', 'var')
    a = lmin / sqrt(3) * 0.5;
    fprintf('Providing default value for a = %f.\n', a)
end

if ~exist('Ec', 'var')
    Ec = MU / (4 * pi) * log(a / 0.1);
    fprintf('Providing default value for Ec = %f.\n', Ec)
end

if ~exist('rann', 'var')
    rann = 0.9 * lmin;
    fprintf('Providing default value for rann = %f.\n', rann)
end

if ~exist('rntol', 'var')
    rntol = 0.5 * rann;
    fprintf('Providing default value for rntol = %f.\n', rntol)
end

if ~exist('rmax', 'var')
    rmax = 0.5 * rann;
    fprintf('Providing default value for rmax = %f.\n', rmax)
end

if ~exist('mindist', 'var')
    mindist = 2 * rann;
    fprintf('Providing default value for rmax = %f.\n', rmax)
end

if ~exist('printfreq', 'var')
    printfreq = 100;
    fprintf('Providing default value for printfreq = %d.\n', printfreq)
end

if ~exist('plotFreq', 'var')
    plotFreq = 100;
    fprintf('Providing default value for plotFreq = %d.\n', plotFreq)
end

if ~exist('simName', 'var')
    simName = append(date,'_',char(timeofday(datetime)));
    fprintf('Providing default value for simName = %s.\n', simName)
end

if ~exist('saveFreq', 'var')
    saveFreq = 100;
    fprintf('Providing default value for saveFreq = %d.\n', saveFreq)
end

if ~exist('viewangle', 'var')
    viewangle = [-35, 15];
    fprintf('Providing default value for viewangle = [%f, %f].\n', viewangle(1), viewangle(2))
end

if ~exist('n_threads', 'var')
    n_threads = 256;
    fprintf('Providing default value for n_threads = %d.\n', n_threads)
end

if ~exist('para_scheme', 'var')
    para_scheme = 1;
    fprintf('Providing default value for para_scheme = %d.\n', para_scheme)
end

if ~exist('u_dot', 'var')
    u_dot = dx / 160E6;
    fprintf('Providing default value for u_dot = %f.\n', u_dot)
end

if ~exist('f_dot', 'var')
    f_dot = dx / 160E6;
    fprintf('Providing default value for f_dot = %d.\n', f_dot)
end

if ~exist('integrator', 'var')
    integrator = 'int_trapezoid';
    fprintf('Providing default value for integrator = %s.\n', integrator)
end

if ~exist('mobility', 'var')
    mobility = 'mobbcc_bb1b';
    fprintf('Providing default value for mobility = %s.\n', mobility)
end

if ~exist('prescribeDofs', 'var')
    prescribeDofs = 'cantilever_bending';
    fprintf('Providing default value for prescribeDofs = %s.\n', prescribeDofs)
end

if ~exist('boundaryConditions', 'var')
    boundaryConditions = 'deformation_control';
    fprintf('Providing default value for boundaryConditions = %s.\n', boundaryConditions)
end

if ~exist('storeBC', 'var')
    storeBC = 'force_displacement_0';
    fprintf('Providing default value for storeBC = %s.\n', storeBC)
end

%% Initialise structures for storage

mods = struct; % Stores function handles of selected modules
flags = struct; % Stores binary flags
FEM = struct; % Stores FEM variables
tract = struct; % Stores traction data
diffBC = struct; % Stores derivatives of BCs
saveBC = struct; % Stores boundary conditions at each time-step
matpara = struct; % Stores simulation and material parameters

%% Store input variables in structures

% mods
mods.integrator = integrator;
mods.mobility = mobility;
mods.prescribeDofs = prescribeDofs;
mods.boundaryConditions = boundaryConditions;
mods.storeBC = storeBC;

% flags
flags.CUDA_flag = CUDA_flag;
flags.a_trac = a_trac;
flags.calculateR_hat = calculateR_hat;
flags.doremesh = doremesh;
flags.docollision = docollision;
flags.doseparation = doseparation;
flags.dovirtmesh = dovirtmesh;

% FEM
FEM.dx = dx; FEM.dy = dy; FEM.dz = dz;
FEM.mx = mx;

% tracts
tract.n_threads = n_threads;
tract.para_scheme = para_scheme;

% diffBC
diffBC.u_dot = u_dot;
diffBC.f_dot = f_dot;

% matpara
matpara.MU = MU;
matpara.NU = NU;
matpara.a = a;
matpara.Ec = Ec;
matpara.rann = rann;
matpara.rntol = rntol;
matpara.rmax = rmax;
matpara.mindist = mindist;
matpara.lmax = lmax;
matpara.lmin = lmin;
matpara.areamax = areamax;
matpara.areamin = areamin;
matpara.mindist = mindist;
