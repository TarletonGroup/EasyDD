%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Francesco Ferroni
% Defects Group / Materials for Fusion and Fission Power Group
% Department of Materials, University of Oxford
% francesco.ferroni@materials.ox.ac.uk
% January 2014
%
% Flags:
% 7 = fixed / no burger's vector conservation (anywhere valid)
% 6 = on free surface (z=0)
% 67 = outside of domain, fixed.
%
% inhull function - checks whether node is inside or outside (convex!)
% domain, defined by vertices (surfaces are triangulated from vertices)
%
% NB, this is overloaded for both rn and rnnew arrays. When specifying
% flags, always use rn(:,end) or rnnew(:,end).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rnnew,linksnew,connectivitynew,linksinconnectnew]=...
    remesh_surf2(rnnew,linksnew,connectivitynew,linksinconnectnew,vertices,P,fn)

% Beginning of surface remeshing for surface nodes. %%%%%%%%%%%%%%%%%%%%
% Flag all nodes outside of medium.
tess = convhulln(vertices); tol=10E-10;
in = inhull(rnnew(:,1:3), vertices, tess, tol); %node with flag 1 is inside domain, 0 is outside domain.
rn_size = size(rnnew,1);
Index=zeros(rn_size,1);

% create the list of nodes and connections
L1=size(rnnew,1);
nodelist=linspace(1,L1,L1)';
[L2,L3]=size(connectivitynew);
conlist=zeros(L2,(L3-1)/2+1); %(L3-1)/2+1 = max/3 M
conlist(:,1)=connectivitynew(:,1);
for i=1:L2
connumb=conlist(i,1);
conlist(i,2:connumb+1)=linspace(1,connumb,connumb);
end

for i=1:rn_size   

    if rnnew(i,end) == 67 || in(i) == 1
        continue;
    end
    
    %estimate surface centroid nearest to point
    [~,index] = min( (P(:,1) - rnnew(i,1)).^2 + (P(:,2) - rnnew(i,2)).^2 + (P(:,3) - rnnew(i,3)).^2 );
    
    %find surface normal associated with centroid
    Index(i)=index;
     
end

for n=1:L1
    n0=nodelist(n); % the node id which is currently considered.
    
    % Find outside nodes connected with inside node with flag 0
    if rnnew(n,end)==0
        

        numNbrs=conlist(n,1); % the number of neighbor nodes
        %rt=zeros(numNbrs,3);
        
        for i=1:numNbrs
            ii=conlist(n,i+1); % connectionid for this connection
            linkid=connectivitynew(n0,2*ii); % the link id
            posinlink=connectivitynew(n0,2*ii+1);
            n1=linksnew(linkid,3-posinlink); % the neighbor node id
            
            if ~isequal(size(in,1),size(rnnew,1))
            in = inhull(rnnew(:,1:3), vertices, tess, tol);
            end
            
            if in(n1)==0%rnnew(n1,end)==67
                rt=rnnew(n1,1:3)-rnnew(n0,1:3); % calculate the length of the link and its tangent line direction
                Lni = norm(rt);
                if Lni>0
                    rt = rt/Lni;
                elseif Lni<=0
                    disp('Zero length segment on the surface during remesh!');
                    pause;
                end
                
                %For segments constituted by a virtual and a real node, find
                %intersection point of segment with surface plane.
                if n1 > size(Index,1)
                    
                %estimate surface centroid nearest to point
                [~,index] = min( (P(:,1) - rnnew(n1,1)).^2 + (P(:,2) - rnnew(n1,2)).^2 + (P(:,3) - rnnew(n1,3)).^2 );
    
                %find surface normal associated with centroid
                Index(n1)=index;
                end
                plane_id = Index(n1);
                if plane_id == 0 || isempty(plane_id) %update id
                    [~,plane_id] = min( (P(:,1) - rnnew(n1,1)).^2 + (P(:,2) - rnnew(n1,2)).^2 + (P(:,3) - rnnew(n1,3)).^2 );
                end
                plane_point = P(plane_id,:);
                plane_normal = fn(plane_id,:);
                sI = dot( plane_normal , (plane_point - rnnew(n0,1:3)) ) / ...
                    dot( plane_normal , (rnnew(n1,1:3) - rnnew(n0,1:3)) );  %sI is a ratio of sclar product M 
                if sI < 0 || sI > 1
                    disp('Split segment does not intersect plane!');
                end
                newsurfNode = rnnew(n0,1:3) + rt*Lni*sI;
                
                % Split only the flag 0 node connected with the flag 6 node.
                % The velocity of flag 6 node is given to the new node to
                % compute new glide planes of new links in splitnode().
                [rn_size1,rn_size2] = size(rnnew);
                if rn_size2 == 4 %Pre-force evaluation, rn nx3, [rn flag]
                    %splitnode accepts only nx6, so add dummy velocities and then
                    %delete them. Avoids re-writing splitnode.
                    dummy_vel = zeros(rn_size1,3);
                    dummy_rn = [rnnew(:,1:3) , dummy_vel , rnnew(:,end)];
                    vntemp = dummy_rn(n1,4:6);
                    posvel = [newsurfNode, vntemp];
                    [dummy_rn,linksnew,connectivitynew,linksinconnectnew]=splitnode(dummy_rn,linksnew,connectivitynew,linksinconnectnew,n0,ii,posvel);
                    dummy_rn(end,end) = 6; %splitnode always gives flag 0, so we have to change it as 6 for surface node!
                    rnnew = dummy_rn(:,[1:3 end]);
                    linksnew(end,6:8)=linksnew(n0,6:8);
                else %Post-force evaluation, rnnew nx6, [rn vel flag]
                    vntemp = rnnew(n1,4:6);
                    posvel = [newsurfNode, vntemp];
                    [rnnew,linksnew,connectivitynew,linksinconnectnew]=splitnode(rnnew,linksnew,connectivitynew,linksinconnectnew,n0,ii,posvel);
                    rnnew(end,end) = 6; %splitnode always gives flag 0, so we have to change it as 6 for surface node!
                    %testing beta
                    linksnew(end,6:8)=linksnew(n0,6:8);
                end
            end
        end
    end
end

for i=1:rn_size
    %if node is already flagged as virtual fixed node, or is internal
    %skip to next node
    if rnnew(i,end) == 67 || in(i) == 1
        continue;
    end
    
    %extend far away
    rnnew = extend(rnnew,linksnew,i,index,fn);
    if any(any(isnan(rnnew)))
        disp('The new extended node is NaN! See Line 46 in remesh_surf.m');
        pause;
    end
    
end



%find surface nodes connected only with virtual nodes and repair them
%accordingly
for i=1:L1 %length(rnnew)  ET updated- check!
    if rnnew(i,end)~=6
        continue;
    end
    test=0;
    for j=1:connectivitynew(i,1)
        if rnnew(linksnew(connectivitynew(i,2*j),logicswap(connectivitynew(i,2*j+1))),end)~=67 %if the surface node is not connected to a virtual node : break M
            break;
        else
            test=test+1;
        end   
    end
    if test==connectivitynew(i,1) %all linked nodes are virtual (67 flag)
        %estimate surface centroid nearest to point
        if Index(i) == 0 || isnan(Index(i)) %update
            [~,index] = min( (P(:,1) - rnnew(i,1)).^2 + (P(:,2) - rnnew(i,2)).^2 + (P(:,3) - rnnew(i,3)).^2 );
        end
        %extend far away
        rnnew = extend(rnnew,linksnew,i,index,fn);
    end
end

%plotting (debug)
%trisurf(tri, Xb(:,1), Xb(:,2), Xb(:,3), 'FaceColor', 'cyan','FaceAlpha', 0.5);
%hold on
%scatter3(rn(:,1),rn(:,2),rn(:,3),'r');
end

function out = logicswap(in)
if in==1
    out=2;
elseif in==2
    out=1;
else
    disp('Incorrect use of logicswap function');
end
end

function rnnew = extend(rnnew,linksnew,rn_id,plane_id,fn)
    %% Option (1) extend using burgers vector
    %Has issues with screw segments
    %extend point to "infinity", or far away, along Burgers vector of segment.
    p_inf=1e5;
    burg = linksnew(linksnew(:,1)==rn_id,3:5);
    %if the node is connected is part to multiple segments, i.e. a junction
    %with segments of differing burgers vector, take the resultant as the
    %vector along which to extend the node.
    burg = sum(burg,1);
    burg = burg/norm(burg);
    %flip sign if angle between b-vec and plane normal is negative to avoid
    %the point being extended back into the material.
    proj_vector = sign(dot(fn(plane_id,1:3),burg))*burg;
    rnnew(rn_id,1:3) = rnnew(rn_id,1:3) + p_inf*proj_vector; %far away...
    rnnew(rn_id,end) = 67;   %flag node as virtual

    %% Option (2) extend using vector orthogonal to line direction and slip plane
    %has issues with screw segments
%     linkID = linksnew(linksnew(:,1)==rn_id,1:2);
%     l_vec = rnnew(linkID(:,2),1:3) - rnnew(linkID(:,1),1:3);
%     l_vec = sum(l_vec,1)/norm(l_vec);
%     
%     slip_plane = linksnew(linksnew(:,1)==rn_id,6:8);
%     slip_plane = sum(slip_plane,1)/norm(slip_plane);
%     
%     ortho_vec = cross(l_vec,slip_plane);
%     
%     if dot(fn(plane_id,1:3),ortho_vec)<0
%         ortho_vec = -ortho_vec;
%     end
%   
%     rnnew(rn_id,1:3) = rnnew(rn_id,1:3) + 10^5*ortho_vec;
%  
%     %fix point at infinity, using flag=67
%     rnnew(rn_id,end) = 67;
    
    %% Option (3) extend based on projection of surface plane normal onto slip plane
    
    %slip plane normal taken from input file!!! must be correct.
    %could calculate as vector orthogonal to b vec and line vec
    %but will be ill-defined for screw segments.
%     slip_plane = linksnew(linksnew(:,1)==rn_id,6:8);
%     if isempty(slip_plane)
%         slip_plane = linksnew(linksnew(:,2)==rn_id,6:8);
%     end
%     slip_plane = sum(slip_plane,1)/norm(slip_plane);
%     
%     %surface normal
%     surf_plane = fn(plane_id,1:3);
%     %surf_plane = surf_plane/norm(surf_plane);
%     
%     proj_matrix = slip_plane' * slip_plane;
%     proj_vec = surf_plane * ( eye(3) - proj_matrix );
%     
%     %Check (1) - make sure projection is pointing outwards of body
%     if sum(surf_plane.*proj_vec)<0 
%          proj_vec = -proj_vec;
%     end
% 
%     %Check (2) - make sure slip plane and surface normal do no coincide.
%     %scenario unlikely to happen if climb drag is very high...
%     normprojvec = norm(proj_vec);
%     if normprojvec < eps
%          disp('Climb?');
%          proj_vec = surf_plane/norm(surf_plane);
%     end
%     
%     proj_vec = proj_vec/normprojvec;
%     
%     rnnew(rn_id,1:3) = rnnew(rn_id,1:3) + 10^5*proj_vec;
%     rnnew(rn_id,end) = 67;
    
    %% Option (4) extend based on projection of surface normal on to slip plane calculated from burgers vector or nodal velocity
%    
%     p_inf=10^5;   %large number used to desegnate pseudo-infinity. Should be updated to releated to simulated volume size
%     vel_vec=rnnew(rn_id,4:6);   %nodal velocity vector
%     vel_vec=vel_vec/norm(vel_vec);   %normalise velocity
%     b_vec=linksnew(linksnew(:,1)==rn_id,3:5);   %find burgers vector of first connected link
%     if isempty(b_vec)
%         b_vec=linksnew(linksnew(:,2)==rn_id,3:5);
%     end
%     b_vec=b_vec/norm(b_vec);   %normalise burgers vector
%     
%     [n,~]=size(linksnew);   %find line direction
%     for i=1:n
%         if linksnew(i,1)==rn_id
%             lin_vec=rnnew(linksnew(i,2),1:3)-rnnew(rn_id,1:3);
%             break
%         elseif linksnew(i,2)==rn_id
%              lin_vec=rnnew(linksnew(i,1),1:3)-rnnew(rn_id,1:3);
%             break
%         end
%     end
%     
%     lin_vec=lin_vec/norm(lin_vec);   %normalise line direction
%     slip_plane=cross(b_vec,lin_vec);   %calculate slip plane normal for non-screw type
%     check=dot(b_vec,lin_vec);   %check if dislocation is screw type
%     check=sqrt(check*check);
%     if 1-check<eps
%     slip_plane=cross(vel_vec,lin_vec);   %calculate slip plane normal if screw type
%     end
%     slip_plane=slip_plane/norm(slip_plane);   %normalise slip plane
%     
%     %surface normal
%     surf_plane = fn(plane_id,1:3);
%     surf_plane = surf_plane/norm(surf_plane);   %normalise surface normal
%     
%     proj_matrix = slip_plane' * slip_plane;   %calculate projection matrix
%     proj_vec = surf_plane * ( eye(3) - proj_matrix );   %calculate extension direction
%     
%     %Check (1) - make sure extension is pointing outwards of body
%     if sum(surf_plane.*proj_vec)<0 
%          proj_vec = -proj_vec;
%     end
%     
%     proj_vec = proj_vec/norm(proj_vec);   %normalise extension vector
%     
%     rnnew(rn_id,1:3) = rnnew(rn_id,1:3) + p_inf*proj_vec;   %extend node to pseudo-infinity in extension direction
%     rnnew(rn_id,end) = 67;   %flag node as virtual
end