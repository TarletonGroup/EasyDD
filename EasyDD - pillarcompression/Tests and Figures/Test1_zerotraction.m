%% test 1 - zero traction at free surface

load ../restart.mat

vector1 = linspace(vertices(1,1),vertices(2,1),100);
vector2 = linspace(vertices(1,2),vertices(3,2),100);

[X,Y] = meshgrid(vector1,vector2);
Z = zeros(100,100);

surf(X,Y,Z);

