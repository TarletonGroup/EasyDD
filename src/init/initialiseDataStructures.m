function [matpara, mods, flags, FEM, tract, diffBC, saveBC] = ...
    initialiseDataStructures(DD)
    %===============================================================%
    % Daniel Hortelano Roig (29/11/2020)
    % daniel.hortelanoroig@materials.ox.ac.uk 
    
    % Organise input data into structures.
    %===============================================================%
    
    %% Initialise structures
    
    mods = struct; % Stores function handles of selected modules
    flags = struct; % Stores binary flags
    FEM = struct; % Stores FEM variables
    tract = struct; % Stores traction data
    diffBC = struct; % Stores derivatives of BCs
    saveBC = struct; % Stores boundary conditions at each time-step
    matpara = struct; % Stores simulation and material parameters

    %% Store input variables in structures
    
    % mods
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
    
    % simpara
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
    
end
