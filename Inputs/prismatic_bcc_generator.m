function [rn, links, normal, b_vec] = prismatic_bcc_generator(n_source, ...
                                   d_source, crys_struct, a_m, dx, dy,  ...
                                   dz, p_lim, l_source)
%=========================================================================%
%-------------------------------------------------------------------------%
% Daniel Celis Garza
% 04/07/2018
%-------------------------------------------------------------------------%
%
% Prismatic bcc loop generator.
%
%=========================================================================%

if strcmp(crys_struct,'bcc')
    disp('Crystal structure recognized');
else
    disp('Crystal structure not recognized. Aborting');
    return;
end

%% BCC Slip systems.
ALLslips = [1 0 1 -1 1 1;
            1 0 1 1 1 -1;
            1 0 -1 1 1 1;
            1 0 -1 1 -1 1;
            1 1 0 -1 1 1;
            1 1 0 1 -1 1;
            1 -1 0 1 1 1;
            1 -1 0 1 1 -1;
            0 1 1 1 -1 1;
            0 1 1 1 1 -1;
            0 1 -1 1 1 1;
            0 1 -1 -1 1 1;];

%% Dislocation with external segments.
i = 1; %to 8
dx = 100;
dy = 50;
dz = 25;

% Slip plane and bugers vector.
slip_n = ALLslips(i,1:3);
b_vec  = ALLslips(i,4:6);

% seg_t1 is the in-plane segment. seg_t2 is the out of plane segment.
seg_t1 = cross(slip_n, b_vec);
seg_t1 = seg_t1/norm(seg_t1);
seg_t2 = slip_n;
seg_t2 = seg_t2/norm(seg_t2);

% Midpoint of the cantilever. Dislocation will pass through it.
mid = [0.5*dx 0.5*dy 0.5*dz];

% Intersect seg_1t with plane.
p_0 = [    0 0.5*dy 0.5*dz;... % midpoint of min(x), yz-plane, face 1
          dx 0.5*dy 0.5*dz;... % midpoint of max(x), yz-plane, face 2
      0.5*dx      0 0.5*dz;... % midpoint of min(y), xz-plane, face 4
      0.5*dx     dy 0.5*dz;... % midpoint of max(y), xz-plane, face 3
      0.5*dx 0.5*dy      0;... % midpoint of min(z), xy-plane, face 5
      0.5*dx 0.5*dy    dz];    % midpoint of max(z), xy-plane, face 6

l_0 = mid;

vertices = [ 0,  0,  0;...
            dx,  0,  0;...
             0, dy,  0;...
            dx, dy,  0;...
             0,  0, dz;...
            dx,  0, dz;...
             0, dy, dz;...
            dx, dy, dz];
        
normals = [-1  0  0;... % min(x), yz-plane, face 1
            1  0  0;... % max(x), yz-plane, face 2
            0 -1  0;... % min(y), xz-plane, face 4
            0  1  0;... % max(y), xz-plane, face 3
            0  0 -1;... % min(z), xy-plane, face 5
            0  0  1];   % max(z), xy-plane, face 6
        
% d = dot(p_0 - l_0, normal)/dot(l, normal)
% intersect = d*l + l_0
l = seg_t1;
int_cntr = 0;
r1r2 = zeros(2,3);
for j = 1:6
    % We only want two intercepts, if it exits at a corner then we pick the
    % exit at the first plane found.
    if (int_cntr > 2)
        break
    end %if
    d = dot((p_0(j,:) - l_0), normals(j,:))/dot(l, normals(j,:));
    intersect = d*l + l_0
    if(any(intersect < 0) || intersect(1) > dx || intersect(2) > dy || intersect(3) > dz)
        continue
    end %if
    int_cntr = int_cntr + 1;
    r1r2(int_cntr,:) = intersect;
end %for

rn = zeros(2,4);
links = zeros(1,8);

rn(1,:) = [r1r2(1,:),7];
rn(2,:) = [r1r2(2,:),7];

links(1,1:2) = [1, 2];
links(1,3:5) = b_vec;
links(1,6:8) = slip_n;

plotnodes(rn,links,0,vertices)
end %function