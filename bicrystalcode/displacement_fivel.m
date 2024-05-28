%Fivel method for calculating dislocations

% A = [Ax,Ay,Az] positions of first node
% B = [Bx,By,Bz] positions of second node
% b = [bx,by,bz] burgers vector of segment
% n = [nx,ny,nz] segment slip plane normal
% l = [lx,ly,lz] line segment direction
% u = [ux,uy,uz] arbitrary vector, origin [0,0,0] (constant for all segs!)
% Cp = [Cpx,Cpy,Cpz] intersection point between slip plane and arbitrary
% vector
% Ap = [Apx,Apy,Apz] projection of Cp onto line collinear with b passing
% through A
% Bp = [Bpx,Bpy,Bpz] projection of Cp onto line collinear with b passing
% through B
    
function [Ux, Uy, Uz] = displacement_fivel(x0,segments,NU)
    
    x0_length = size(x0,1);
    utilda_array = zeros(x0_length,3);

    for i=1:size(segments,1)
        b = segments(i,3:5);
        A = segments(i,6:8);
        B = segments(i,9:11);
        l = B-A;
%         u = [1, 0, 0]; % u.n must be non zero for all slip planes normals n
        u = [sqrt(3),0.5,2+sqrt(3)]; % is this the arbitary vector
        u_origin = [0,0,0];

        % Calculate slip plane of segment
        n = cross(l,b);
        normn = norm(n);
        n = n./normn; %unit vector        
      
        if any(isnan(n)) || all(n==0) || normn < 1E-10
            %screw segment - slip plane of screw is indetermined
            %assign n equivalent to that given in links
            n = segments(i,12:14);
            %Looks at adjacent segment. Breaks if two consecutive screws are present!
            % connects=find(segments(i,1)==segments(:,2));
            % b_conn=segments(connects,3:5);
            % l_conn=segments(connects,9:11)-segments(connects,6:8);
            % n=cross(b_conn,l_conn,2);
            % for j=1:size(n,1)
            %    n(j,:)=n(j,:)./norm(n(j,:));
            % end
            % logic=not(any(isnan(n)));
            % n=n(logic,:);
            
             n = segments(i,12:14); % lets use n from links.
             normn =norm(n);
             n = n./normn;             
             if any(isnan(n)) || all(n==0) || normn < 1E-10
                 disp('Incorrect input n in links');
                 pause;
             end
        end
        
        % the sign of n does not matter for finding Cp :)
      
        % Intersection between arbitrary vector u (origin 0) and plane
        [Cp,check] = plane_line_intersect(n,A,u_origin,u);

        if check == 2 || check == 0
            disp('Choose different arbitrary vector for displacement calculations');
        end

        % Find orthogonal projection of closure point Cp on the line collinear with
        % b which passes through a node.
        Ap = A + (sum((Cp-A).*b)/sum(b.*b))*b;
        Bp = B + (sum((Cp-B).*b)/sum(b.*b))*b;
        
        % Calculate displacement due to the two subtriangles associated with each
        % segment
        %
        %   A --> B
        %   ^ \   |
        %   |  \  |
        %   |   \ | 
        %   Ap<--Bp      . point where displacement is needed
        %
        % In other words, of triangles Bp->A->B and Ap->A->Bp with point of
        % interest

        for j=1:x0_length
            coord = x0(j,:);
            
%            utilda = displacementmex(coord,Bp,A,B,b,n,NU) + displacementmex(coord,Ap,A,Bp,b,n,NU);
%              utilda1 = displacement_et(coord,Bp,A,B,b,NU) + displacement_et(coord,Ap,A,Bp,b,NU);
             utilda = displacementmex_et(coord,Bp,A,B,b,NU) + displacementmex_et(coord,Ap,A,Bp,b,NU);
            % == displacementmex(coord,A,B,Ap,b,n,NU) + displacementmex(coord,Ap,B,Bp,b,n,NU)
            if any(not(isfinite(utilda)))
%                 tol=10e-5;
%                 displ=10e-2;
%                 if sum(abs(abs(Ap)-abs(Bp))) < tol
                    %In this case, the A' and B' points are effectively the same,
                    %meaning that it's a pure screw segment. 
                    utilda = [0; 0; 0]; %pure screw has zero displacement
%                 else
%                     %Include check of collinearity between "coord" and AB,BBp,BpAp,AAp,ABp
%                     %as this may produce numerical issues, like +/-inf or NaN.
%                     % -> Add a small distance to coord, away from line!
%                     btol = norm(b)*10E-4;
%                     %Checks if three points are collinear.
%                     AAp_logic = rank([coord-Ap;A-Ap],btol) < 2;
%                     %AB_logic = rank([coord-B;A-B],btol) < 2;
%                     ABp_logic = rank([coord-Bp;A-Bp],btol) < 2;
%                     ApBp_logic = rank([coord-Bp;Ap-Bp],btol) < 2;
%                     BBp_logic = rank([coord-Bp;B-Bp],btol) < 2;
%                     disp('Singular utildas. Correcting...');
%                     %In this case, the segment is not pure screw, so 2 Barnett
%                     %triangles need to be evaluated to correctly calculate its utilda
%                     %contribution.
%                     if ABp_logic == 1 || BBp_logic == 1 %|| AB_logic == 1
%                         %In this case, the coordinate is collinear to some segments of
%                         %the Barnett triangle B'AB. Move the evaluation point slightly
%                         %away to avoid numerical issues (Inf/NaN).
%                         utilda = displacementmex(coord,Ap,A,Bp,b,n,NU) + displacementmex(coord+displ*n,Bp,A,B,b,n,NU);
%                         %utilda = displacementmex(coord,Ap,A,Bp,b,n,NU);
%                     elseif AAp_logic == 1 || ApBp_logic == 1  %|| ABp_logic == 1
%                         %In this case, the coordinate is collinear to some segments of
%                         %the Barnett triangle B'AA'. Move the evaluation point slightly
%                         %away to avoid numerical issues (Inf/NaN).
%                         utilda = displacementmex(coord,Bp,A,B,b,n,NU) + displacementmex(coord+displ*n,Ap,A,Bp,b,n,NU);
%                         %utilda = displacementmex(coord,Bp,A,B,b,n,NU);
%                     else
%                        disp('Corrections failed. Zeroed for now...');
%                        utilda = [0; 0; 0];
%                 %    end
%                 end
            end
        
            utilda_array(j,1:3) = utilda_array(j,1:3) + utilda';
        end
    end
    
    Ux = utilda_array(:,1);
    Uy = utilda_array(:,2);
    Uz = utilda_array(:,3);
end

function [I,check]=plane_line_intersect(n,V0,P0,P1)
%plane_line_intersect computes the intersection of a plane and a segment(or
%a straight line)
% Inputs: 
%       n: normal vector of the Plane 
%       V0: any point that belongs to the Plane 
%       P0: end point 1 of the segment P0P1
%       P1:  end point 2 of the segment P0P1
%
% Outputs:
%      I    is the point of interection 
%     Check is an indicator:
%      0 => disjoint (no intersection)
%      1 => the plane intersects P0P1 in the unique point I
%      2 => the segment lies in the plane
%      3=>the intersection lies outside the segment P0P1

I=[0 0 0];
u = P1-P0;
%w = P0 - V0;
w = V0 - P0;
D = sum(n.*u); %dot product
%N = -sum(n.*w); %dot product
N = sum(n.*w); %dot product
check=0;
if abs(D) < 10^-7        % The segment is parallel to plane
        if N == 0           % The segment lies in plane
            check=2;
            return
        else
            check=0;       %no intersection
            return
        end
end

%compute the intersection parameter
sI = N / D;
I = P0+ sI.*u;

if (sI < 0 || sI > 1)
    check= 3;          %The intersection point  lies outside the segment, so there is no intersection
else
    check=1;
end

end