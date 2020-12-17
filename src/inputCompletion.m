%===============================================================%
% Ensure we have all input variables.
% Daniel Celis Garza, Aug 2020
%===============================================================%
if ~exist('amag', 'var')
    amag = 3.18e-4;
    fprintf('Providing default value for amag = %f.\n', amag)
end

if ~exist('mumag', 'var')
    mumag = 1.6e5; % MPa only for plotting.
    fprintf('Providing default value for mumag = %f.\n', mumag)
end

if ~exist('CRYSTAL_STRUCTURE', 'var')
    CRYSTAL_STRUCTURE = 'bcc';
    fprintf('Providing default value for CRYSTAL_STRUCTURE = %s.\n', CRYSTAL_STRUCTURE)
end

if ~exist('NUM_SOURCES', 'var')
    NUM_SOURCES = 1;
    fprintf('Providing default value for NUM_SOURCES = %f.\n', NUM_SOURCES)
end

if ~exist('DIST_SOURCE', 'var')
    DIST_SOURCE = 0.5 / amag;
    fprintf('Providing default value for DIST_SOURCE = %f.\n', DIST_SOURCE)
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
    my = round(mx * dy / dx); % # elements in y direction
    my = max(my, 1);
    mz = round(mx * dz / dx); % elements in z direction
    mz = max(mz, 1);
    fprintf('Providing default value for mx = %f.\n', mx)
end

if ~exist('loading', 'var')
    loading = @displacementControl;
    fprintf('Providing default value for loading:\n'),
    functions(loading)
end

if ~exist('calculateLoading', 'var')
    calculateLoading = @constantLoading;
    fprintf('Providing default value for calculateLoading:\n')
    functions(calculateLoading)
end

if ~exist('vertices', 'var')
    vertices = [0, 0, 0; ...
                dx, 0, 0; ...
                0, dy, 0; ...
                dx, dy, 0; ...
                0, 0, dz; ...
                dx, 0, dz; ...
                0, dy, dz; ...
                dx, dy, dz];
    fprintf('Providing default value for vertices.\n')
end

if ~exist('faces', 'var')
    faces = [1, 3, 4, 2; % Faces of cuboid as defined by vertices
        5, 6, 8, 7;
        2, 4, 8, 6;
        1, 5, 7, 3;
        1, 2, 6, 5;
        3, 7, 8, 4];
    fprintf('Providing default value for faces.\n')
end

if ~exist('MU', 'var')
    MU = 1;
    fprintf('Providing default value for MU = %f.\n', MU)
end

if ~exist('NU', 'var')
    NU = 0.305;
    fprintf('Providing default value for NU = %f.\n', NU)
end

if ~exist('mobility', 'var')
    mobility = @mobbcc_bb1b;
    fprintf('Providing default value for mobility:\n')
    functions(mobility)
end

if ~exist('Bcoeff', 'var')
    Bcoeff = struct('screw', 10, 'edge', 1, 'climb', 1e10, 'line', 1e-4);
    fprintf('Providing default value for Bcoeff.screw = %f.\n', Bcoeff.screw)
    fprintf('Providing default value for Bcoeff.edge = %f.\n', Bcoeff.edge)
    fprintf('Providing default value for Bcoeff.climb = %f.\n', Bcoeff.climb)
    fprintf('Providing default value for Bcoeff.line = %f.\n', Bcoeff.line)
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

if ~exist('integrator', 'var')
    integrator = @int_trapezoid;
    fprintf('Providing default value for integrator:\n')
    functions(integrator)
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

if ~exist('printfreq', 'var')
    printfreq = 100;
    fprintf('Providing default value for printfreq = %d.\n', printfreq)
end

if ~exist('plotFreq', 'var')
    plotFreq = 100;
    fprintf('Providing default value for plotFreq = %d.\n', plotFreq)
end

if ~exist('plotFlags', 'var')
    plotFlags = struct('nodes', true, 'secondary', false);
    fprintf('Providing default value for plotFlags.nodes = %d.\n', plotFlags.nodes)
    fprintf('Providing default value for plotFlags.secondary = %d.\n', plotFlags.secondary)
end

if ~exist('simName', 'var')
    simName = date;
    fprintf('Providing default value for simName = %s.\n', simName)
end

if ~exist('saveFreq', 'var')
    saveFreq = 100;
    fprintf('Providing default value for saveFreq = %d.\n', saveFreq)
end

if ~exist('plim', 'var')
    plim = max([dx, dy, dz]) / amag;
    fprintf('Providing default value for plim = %f.\n', plim)
end

if ~exist('viewangle', 'var')
    viewangle = [-35, 15];
    fprintf('Providing default value for viewangle = [%f, %f].\n', viewangle(1), viewangle(2))
end

if ~exist('calculateTractions', 'var')
    calculateTractions = @calculateAnalyticTractions;
    fprintf('Providing default value for calculateTractions:\n')
    functions(calculateTractions)
end

if ~exist('CUDA_flag', 'var')
    CUDA_flag = false;
    fprintf('Providing default value for CUDA_flag = %d.\n', CUDA_flag)
end

if ~exist('n_threads', 'var')
    n_threads = 256;
    fprintf('Providing default value for n_threads = %d.\n', n_threads)
end

if ~exist('para_scheme', 'var')
    para_scheme = 1;
    fprintf('Providing default value for para_scheme = %d.\n', para_scheme)
end

if ~exist('sign_u_dot', 'var')
    sign_u_dot = -1;
    fprintf('Providing default value for sign_u_dot= %d.\n', sign_u_dot)
end

if ~exist('sign_f_dot', 'var')
    sign_f_dot = -1;
    fprintf('Providing default value for sign_f_dot= %d.\n', sign_f_dot)
end

if ~exist('u_dot', 'var')
    u_dot = dx / 160E6;
    fprintf('Providing default value for u_dot = %f.\n', u_dot)
end

if ~exist('u_bar', 'var')
    u_bar = 0;
    fprintf('Providing default value for u_bar = %f.\n', u_bar)
end

if ~exist('f_dot', 'var')
    f_dot = dx / 160E6;
    fprintf('Providing default value for f_dot = %d.\n', f_dot)
end

if ~exist('f_bar', 'var')
    f_bar = 0;
    fprintf('Providing default value for f_bar = %f.\n', f_bar)
end

if ~exist('f_hat', 'var')
    f_hat = 0;
    fprintf('Providing default value for f_hat = %f.\n', f_hat)
end

if ~exist('f_tilda', 'var')
    f_tilda = 0;
    fprintf('Providing default value for f_tilda = %f.\n', f_tilda)
end

if ~exist('u_hat', 'var')
    u_hat = 0;
    fprintf('Providing default value for u_hat = %f.\n', u_hat)
end

if ~exist('u_tilda', 'var')
    u_tilda = 0;
    fprintf('Providing default value for u_tilda = %f.\n', u_tilda)
end

if ~exist('r_hat', 'var')
    r_hat = 0;
    fprintf('Providing default value for r_hat = %f.\n', r_hat)
end

if ~exist('calculateLoadingFunctionArgs', 'var')
    calculateLoadingFunctionArgs = 0;
    fprintf('Providing default value for calculateLoadingFunctionArgs = %f\n.', calculateLoadingFunctionArgs)
end

if ~exist('simType', 'var')
    simType = @cantileverBending;
    fprintf('Providing default value for simType:\n')
    functions(simType)
end

if ~exist('Fsim', 'var')
    Fsim = zeros(1e6, 1);
    fprintf('Initialising Fsim as 1e6 zeros.\n')
end

if ~exist('Usim', 'var')
    Usim = zeros(1e6, 1);
    fprintf('Initialising Usim as 1e6 zeros.\n')
end

if ~exist('t', 'var')
    t = zeros(1e6, 1);
    fprintf('Initialising t as 1e6 zeros.\n')
end

if ~exist('dt', 'var')
    dt = dt0;
    fprintf('Providing default value for dt = %f.\n', dt)
end

if ~exist('rotMatrix', 'var')
    rotMatrix = [];
    fprintf('Using standard coordinates, x = [1 0 0], y = [0 1 0], z = [0 0 1]\n')
end
