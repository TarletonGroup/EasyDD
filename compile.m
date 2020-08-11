function CUDA_flag = compile(CUDA_flag)
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
name = sprintf("SegSegForcesMex.%s", ext);
file = dir(name);
if ~isfile(name) || ~isfile(file.name) && days(file.date - datetime('now')) > 30
    mex -v COPTIMFLAGS="-o3 -oy -use_fast_math -DNDEBUG" SegSegForcesMex.c
end

% Stress tue to segments (field point stress).
name = sprintf("StressDueToSegsMex.%s", ext);
file = dir(name);
if ~isfile(name) || ~isfile(file.name) && days(file.date - datetime('now')) > 30
    mex -v COPTIMFLAGS="-o3 -oy -use_fast_math -DNDEBUG" StressDueToSegsMex.c
end

% Stress tue to segments (field point stress).
file = dir(sprintf("MinDistCalcMex.%s", ext));
if ~isfile("MinDistCalcMex.c") || ~isfile(file.name) && days(file.date - datetime('now')) > 30
    mex -v COPTIMFLAGS="-o3 -oy -use_fast_math -DNDEBUG" MinDistCalcMex.c
end

% Collision checker.
name = sprintf("CollisionCheckerMex.%s", ext);
file = dir(name);
if ~isfile(name) || ~isfile(file.name) && days(file.date - datetime('now')) > 30
    mex -v COPTIMFLAGS="-o3 -oy -use_fast_math -DNDEBUG" CollisionCheckerMex.c
end

% Analytic surface tractions.
name = sprintf("NodalSurfaceForceLinearRectangleMex.%s", ext);
file = dir(name);
if ~isfile(name) || ~isfile(file.name) && days(file.date - datetime('now')) > 30
    mex -v COPTIMFLAGS="-o3 -oy -use_fast_math -DNDEBUG" NodalSurfaceForceLinearRectangleMex.c
end


if CUDA_flag
    try
        % Analytic surface tractions CUDA.
        name = sprintf("NodalSurfaceForceLinearRectangleMexCUDA.%s", ext);
        file = dir(name);
        if ~isfile(name) || ~isfile(file.name) && days(file.date - datetime('now')) > 30
            mexcuda -v COMPFLAGS="-Xptxas" COPTIMFLAGS="-o3 -oy -use_fast_math -DNDEBUG" NodalSurfaceForceLinearRectangleMexCUDA.cu
        end
        
        % SegSeg forces CUDA.
        name = sprintf("SegForceNBodyCUDADoublePrecision.%s", ext);
        file = dir(name);
        if ~isfile(name) || ~isfile(file.name) && days(file.date - datetime('now')) > 30
            system('nvcc -ptx -o3 -oy -use_fast_math -DNDEBUG SegForceNBodyCUDADoublePrecision.cu');
        end     
    catch err
        CUDA_flag = false;
        disp("No CUDA compiler found. Using serial implementation. CUDA_flag set to false.")
    end
end

end