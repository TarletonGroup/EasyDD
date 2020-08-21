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
    dy = 5 / amag; % 5 micron
    fprintf('Providing default value for dz = %f.\n', dz)
end

if ~exist('mx', 'var')
    mx = 20;
    fprintf('Providing default value for mx = %f.\n', mx)
end

if ~exist('loading', 'var')
    loading = 1;
    fprintf('Providing default value for loading = %d.', loading)
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
    mobility = 'mobbcc_bb1b';
    fprintf('Providing default value for mobility = %s.\n', mobility)
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
    fprintf('Providing default value for doremesh = %f.\n', doremesh)
end

if ~exist('docollision', 'var')
    docollision = 1; %flat set to 0 or 1 that turns collision detection off or on
    fprintf('Providing default value for docollision = %f.\n', docollision)
end

if ~exist('doseparation', 'var')
    doseparation = 1; %flat set to 0 or 1 that turns splitting algorithm for highly connected node off or on
    fprintf('Providing default value for doseparation = %f.\n', doseparation)
end

if ~exist('dovirtmesh', 'var')
    dovirtmesh = 1; %flat set to 0 or 1 that turns remeshing of virtual nodes off or on
    fprintf('Providing default value for dovirtmesh = %f.\n', dovirtmesh)
end

if ~exist('dt0', 'var')
    dt0 = 1e5;
    fprintf('Providing default value for dt0 = %f.\n', dt0)
end

if ~exist('', 'var')
    fprintf('Providing default value for = %f.\n', )
end

if ~exist('', 'var')
    fprintf('Providing default value for = %f.\n', )
end

if ~exist('', 'var')
    fprintf('Providing default value for = %f.\n', )
end

if ~exist('', 'var')
    fprintf('Providing default value for = %f.\n', )
end

if ~exist('', 'var')
    fprintf('Providing default value for = %f.\n', )
end

if ~exist('', 'var')
    fprintf('Providing default value for = %f.\n', )
end

if ~exist('', 'var')
    fprintf('Providing default value for = %f.\n', )
end

if ~exist('', 'var')
    fprintf('Providing default value for = %f.\n', )
end

if ~exist('CUDA_flag', 'var')
    CUDA_flag = false;
end
