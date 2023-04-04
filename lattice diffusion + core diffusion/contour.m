clear
clc
S = load('./output/1.dat');
x1 = S(:,1);
y1 = S(:,2);
z1 = S(:,3);
x1 = x1';
y1 = y1';
z1 = z1';

S = load('./output/1000.dat');
x2 = S(:,1);
y2 = S(:,2);
z2 = S(:,3);
x2 = x2';
y2 = y2';
z2 = z2';

S = load('./output/2000.dat');
x3 = S(:,1);
y3 = S(:,2);
z3 = S(:,3);
x3 = x3';
y3 = y3';
z3 = z3';

S = load('./output/3000.dat');
x4 = S(:,1);
y4 = S(:,2);
z4 = S(:,3);
x4 = x4';
y4 = y4';
z4 = z4';

S = load('./output/4000.dat');
x5 = S(:,1);
y5 = S(:,2);
z5 = S(:,3);
x5 = x5';
y5 = y5';
z5 = z5';


xi = 0 : 5: 625;
yi = -320 : 5 : 320;

[XI, YI] = meshgrid(xi, yi);

% subplot(2,2,1)
% Z1 = griddata(x, y, z, X1, Y1,'linear');
% mesh(X1, Y1, Z1)
% surf(X1, Y1, Z1)
% shading interp;
% title('linear')
% 
% subplot(2,2,2)
% Z1 = griddata(x, y, z, X1, Y1,'nearest');
% mesh(X1, Y1, Z1)
% surf(X1, Y1, Z1)
% shading interp;
% title('nearest')
% 
% subplot(2,2,3)
% Z1 = griddata(x, y, z, X1, Y1,'cubic');
% mesh(X1, Y1, Z1)
% surf(X1, Y1, Z1)
% shading interp;
% title('cubic')
% 
% subplot(2,2,4)
Z1 = griddata(x1, y1, z1, XI, YI,'v4');
Z2 = griddata(x2, y2, z2, XI, YI,'v4');
Z3 = griddata(x3, y3, z3, XI, YI,'v4');
Z4 = griddata(x4, y4, z4, XI, YI,'v4');
Z5 = griddata(x5, y5, z5, XI, YI,'v4');

% mesh(XI, YI, Z1)
surf(XI, YI, Z1)
hold on

% mesh(XI, YI, Z2)
% surf(XI, YI, Z2)
hold on

% mesh(XI, YI, Z3)
surf(XI, YI, Z3)
hold on

surf(XI, YI, Z4)
hold on

surf(XI, YI, Z5)
% 
shading interp;
hold off

% title('v4')