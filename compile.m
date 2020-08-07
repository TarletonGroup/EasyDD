function compile(CUDA_flag)
%%=======================================================================%%
%-------------------------------------------------------------------------%
% Compile MEX and MEXCUDA files.
% Compiles files if they have not been compiled of if they were compiled
% over 30 days ago.
%-------------------------------------------------------------------------%
%=========================================================================%
% Inputs
%=========================================================================%
% CUDA_flag := flag in case CUDA codes required. If 1 compile, else do not
% compile.
%=========================================================================%
% Dummy variables.
%=========================================================================%
% ext := compiled file extension, is dependent on OS.
%-------------------------------------------------------------------------%
% If any new C or CUDA files are added to EasyDD, place them here.
%%=======================================================================%%
if ~exist('CUDA_flag','var')
    CUDA_flag = false;
end

ext = mexext;

% Seg seg forces.
file = dir(sprintf("SegSegForcesMex.%s", ext));
if ~isfile("SegSegForcesMex.c") || ~isfile(file.name) && days(file.date - datetime('now')) > 30
    mex -v COPTIMFLAGS="-o3 -oy -DNDEBUG" SegSegForcesMex.c   
end

% Stress tue to segments (field point stress).
file = dir(sprintf("StressDueToSegsMex.%s", ext));
if ~isfile("StressDueToSegsMex.c") || ~isfile(file.name) && days(file.date - datetime('now')) > 30
    mex -v COPTIMFLAGS="-o3 -oy -DNDEBUG" StressDueToSegsMex.c   
end

% Stress tue to segments (field point stress).
file = dir(sprintf("MinDistCalcMex.%s", ext));
if ~isfile("MinDistCalcMex.c") || ~isfile(file.name) && days(file.date - datetime('now')) > 30
    mex -v COPTIMFLAGS="-o3 -oy -DNDEBUG" MinDistCalcMex.c   
end


% mex -v COPTIMFLAGS="-O3 -Oy -DNDEBUG" MinDistCalcMex.c
% mex -v COPTIMFLAGS="-O3 -Oy -DNDEBUG" CollisionCheckerMex.c
% mex -v COPTIMFLAGS="-O3 -Oy -DNDEBUG" mobbcc1mex.c
% mex -v COPTIMFLAGS="-O3 -Oy -DNDEBUG" displacementmex_et.c
% mex -v COPTIMFLAGS="-O3 -Oy -DNDEBUG" CreateInputMex.c %CollisionMarielle
% mex -v COPTIMFLAGS="-O3 -Oy -DNDEBUG" CollisionCheckerMexMariellebis.c %CollisionMarielle
% mex -v COPTIMFLAGS="-O3 -Oy -DNDEBUG" ./nodalForce/nodal_surface_force_linear_rectangle_mex.c
% mexcuda -v COMPFLAGS="-Xptxas -v -arch=compute_60 -code=sm_60 -use_fast_math" ./nodalForce/nodal_surface_force_linear_rectangle_mex_cuda.cu
end