function [rn, links, slip_n, b_vec] = prismatic_bcc_generator(slips, ...
                                                              dx, dy, dz, ...
                                                              vertices, fem_n)
%=========================================================================%
%-------------------------------------------------------------------------%
% Daniel Celis Garza
% 04/07/2018
%-------------------------------------------------------------------------%
%
% Prismatic bcc loop generator.
%
%=========================================================================%

%% FEM parameters.
ALLslips = [1  0  1 -1  1  1;
            1  0  1  1  1 -1;
            1  0 -1  1  1  1;
            1  0 -1  1 -1  1;
            1  1  0 -1  1  1;
            1  1  0  1 -1  1;
            1 -1  0  1  1  1;
            1 -1  0  1  1 -1;
            0  1  1  1 -1  1;
            0  1  1  1  1 -1;
            0  1 -1  1  1  1;
            0  1 -1 -1  1  1];
        
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
        
%% Dislocations.
n_slips = size(slips, 1);
r1r2 = zeros(2, 3);
rn = zeros(2*n_slips, 4);
links = zeros(1*n_slips, 8);

for i = 1: n_slips
    % Slip plane and bugers vector.
    slip_n = ALLslips(slips(i),1:3);
    b_vec  = ALLslips(slips(i),4:6);
    b_vec = sqrt(3)/2*b_vec/norm(b_vec);

    % seg_t1 is the in-plane segment. seg_t2 is the out of plane segment.
    dln = cross(slip_n, b_vec);
    dln = dln/norm(dln);

    % d = dot(p_0 - l_0, normal)/dot(l, normal)
    % intersect = d*l + l_0
    l = dln;
    int_cntr = 0;

    for j = 1:6
        % We only want two intercepts, if it exits at a corner then we pick the
        % exit at the first plane found.
        if (int_cntr > 2)
            break
        end %if
        d = dot((p_0(j,:) - l_0), fem_n(j,:))/dot(l, fem_n(j,:));
        intersect = d*l + l_0;
        if(any(intersect < 0) || intersect(1) > dx || intersect(2) > dy || intersect(3) > dz)
            continue
        end %if
        int_cntr = int_cntr + 1;
        r1r2(int_cntr,:) = intersect;
    end %for
    
    rn((i-1)*2+1, :) = [r1r2(1, :) 7]; 
    rn((i-1)*2+2, :) = [r1r2(2, :) 7];
    links((i-1)*1+1, 1:2) = [(i-1)*1+1, (i-1)*1+2];
    
    links(((i-1)*1+1), 3:5) = repmat(b_vec , 1, 1);
    links(((i-1)*1+1), 6:8) = repmat(slip_n, 1, 1);
end %for

end %function