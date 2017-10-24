% Parallel dislocation segment to plane.
x1 = [0.5 0. 0.5];
x2 = [0.5 1.0 1.5];
%x3 = [0.0 0.0 0.0];
%x4 = [1.0 0.0 0.0];
%x5 = [0.0 1.0 1.0];
%x6 = [1.0 1.0 1.0];
%x1 = [0.5 0.0 0.75];
%x2 = [0.25 1.2 1.5];
x3 = [0.0 0.0 0.0];
x4 = [1.0 0.0 0.0];
x5 = [0.0 1.0 1.0];
x6 = [1.0 1.0 1.0];
b = [-3. 5. 4.]/10.;
mu = 0.8;
nu = 0.5;
a = 0.01;
time = 0.0;
time2 = 0.0;
for i=1:100
%x1 = rand(1,3);
%x2 = rand(1,3);
tic;
%fprintf('Running MATLAB version...\n');
[fx3,fx4,fx5,fx6,ftot] = NodalSurfForceLinearRectangle2(x1,x2,x3,x4,x5,x6,b,mu,nu,a);
time = time + toc;
ftot
tic;
%fprintf('Running C version (Daniel)...\n');
[fx3mex,fx4mex,fx5mex,fx6mex,ftotmex] = nodal_surface_force_linear_rectangle(x1,x2,x3,x4,x5,x6,b,mu,nu,a);
time2 = time2 + toc;
ftotmex'
%fprintf('accuracy = [%f, %f, %f]\n', ftot(1)/ftotmex(1), ftot(2)/ftotmex(2), ftot(3)/ftotmex(3));
end
fprintf('speed up = %3.6f\n', time/time2);