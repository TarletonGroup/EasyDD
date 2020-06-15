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
    remesh_surf(rnnew,linksnew,connectivitynew,linksinconnectnew,vertices,P,fn)

% Beginning of surface remeshing for surface nodes. %%%%%%%%%%%%%%%%%%%%
% Flag all nodes outside of medium.
tess = convhulln(vertices); tol=10E-10;
in = inhull(rnnew(:,1:3), vertices, tess, tol); %node with flag 1 is inside domain, 0 is outside domain.
rn_size = size(rnnew,1);
Index=zeros(rn_size,1);

% Create the list of nodes and connections
L1=size(rnnew,1);
nodelist=linspace(1,L1,L1)';
[L2,L3]=size(connectivitynew);
conlist=zeros(L2,(L3-1)/2+1); %(L3-1)/2+1 = max/3 M
conlist(:,1)=connectivitynew(:,1);
for i=1:L2
connumb=conlist(i,1);
conlist(i,2:connumb+1)=linspace(1,connumb,connumb);
end

%Find newly exited nodes and flag them as such
for n=1:L1
    if rnnew(n,end)==0 && in(n)==0
        rnnew(n,end)=65;
    end
end

%Halt dislocations attempting to exit the fixed end (assuming fixed end is the y/z plane at x=0)
for n=1:L1
   if rnnew(n,end)==65 && rnnew(n,1)<=0
       vec=[0,0,0];
       connumb=connectivitynew(n,1);
       for m=1:connumb
           vec=vec+rnnew(linksnew(connectivitynew(n,2*m),3-connectivitynew(n,2*m+1)),1:3)-rnnew(n,1:3);
       end
       vec=rnnew(n,1).*(vec/vec(1,1));
       rnnew(n,1:3)=rnnew(n,1:3)-vec;
       rnnew(n,end)=7;
       if any(isnan(rnnew(n,1:3)))
           pause;
       end
   end
end

%Create surface nodes for newly exited nodes
for n=1:L1
    n0=nodelist(n); % the node id which is currently considered.

    % Find outside nodes connected with inside node with flag 0
    if rnnew(n0,end)==0

        numNbrs=conlist(n,1); % the number of neighbor nodes

        for i=1:numNbrs
            ii=conlist(n,i+1); % connectionid for this connection
            linkid=connectivitynew(n0,2*ii); % the link id
            posinlink=connectivitynew(n0,2*ii+1);
            n1=linksnew(linkid,3-posinlink); % the neighbor node id

            if rnnew(n1,end)==65
                [rnnew,linksnew,connectivitynew,linksinconnectnew] = gensurfnode2(Index,fn,P,rnnew,linksnew,connectivitynew,linksinconnectnew,n0,n1,i,ii,vertices);
            end
        end
    end
end

%Extend newly exited nodes and flag as virtual
rn_size = size(rnnew,1);
for i=1:rn_size
    %if node is already flagged as virtual fixed node, or is internal
    %skip to next node
    if rnnew(i,end) ~= 65
        continue;
    end

    %estimate surface centroid nearest to point
    [~,index] = min( (P(:,1) - rnnew(i,1)).^2 + (P(:,2) - rnnew(i,2)).^2 + (P(:,3) - rnnew(i,3)).^2 );

    %find surface normal associated with centroid
    Index(i)=index;

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
        %estimate surface centroid nearest to point
        if Index(i) == 0 || isnan(Index(i)) %update
            [~,index] = min( (P(:,1) - rnnew(i,1)).^2 + (P(:,2) - rnnew(i,2)).^2 + (P(:,3) - rnnew(i,3)).^2 );
            Index(i)=index;
        end
        if rnnew(linksnew(connectivitynew(i,2*j),logicswap(connectivitynew(i,2*j+1))),end)==0 %if the surface node is not connected to a virtual node : break M
            break;
        else
            test=test+1;
        end
    end
    if test==connectivitynew(i,1) %all linked nodes are virtual (67 flag)
        %extend far away
        rnnew = extend(rnnew,linksnew,i,index,fn);
    end
end

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
    % Extend along surface normal

    %N.B. This method should not be used with Fivel's displacement code

    p_inf=1e7;
    proj_vector=fn(plane_id,1:3);
    rnnew(rn_id,1:3) = rnnew(rn_id,1:3) + p_inf*proj_vector; %far away...
    rnnew(rn_id,end) = 67;   %flag node as virtual
end

function [rnnew,linksnew,connectivitynew,linksinconnectnew] = gensurfnode(Index,fn,P,rnnew,linksnew,connectivitynew,linksinconnectnew,n0,n1,i,ii)
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
                plane_id = 0;
                if size(Index,1)>=n1
                plane_id = Index(n1);
                end
                if plane_id == 0 || isempty(plane_id) %update id
                    [~,plane_id] = min( (P(:,1) - rnnew(i,1)).^2 + (P(:,2) - rnnew(i,2)).^2 + (P(:,3) - rnnew(i,3)).^2 );
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
                    linksnew(end,6:8)=linksnew(connectivitynew(n0,2),6:8);
                else %Post-force evaluation, rnnew nx6, [rn vel flag]
                    vntemp = rnnew(n1,4:6);
                    posvel = [newsurfNode, vntemp];
                    [rnnew,linksnew,connectivitynew,linksinconnectnew]=splitnode(rnnew,linksnew,connectivitynew,linksinconnectnew,n0,ii,posvel);
                    rnnew(end,end) = 6; %splitnode always gives flag 0, so we have to change it as 6 for surface node!
                    %testing beta
                    linksnew(end,6:8)=linksnew(n0,6:8);
                end

end
function [rnnew,linksnew,connectivitynew,linksinconnectnew] = gensurfnode2(Index,fn,P,rnnew,linksnew,connectivitynew,linksinconnectnew,n0,n1,i,ii,vertices)
   %HY201907030: the sequence of nodes seems wrong

%    faces = [1,2,3,4;                                             %Faces of cuboid as defined by vertices
%             1,2,5,6;
%             1,3,5,7;
%             2,4,6,8;
%             3,4,7,8;
%             5,6,7,8];

%HY20190730: the sequence of nodes correted by HY

   faces = [1,2,4,3;                                             %Faces of cuboid as defined by vertices
            1,2,6,5;
            1,3,7,5;
            2,4,8,6;
            3,4,8,7;
            5,6,8,7];
   rt = rnnew(n1,1:3)-rnnew(n0,1:3);
   line = [rnnew(n1,1:3),rt];
   surfpts = intersectLineMesh3d(line, vertices, faces);
   [~,point_id] = min((surfpts(:,1)-rnnew(n0,1)).^2+(surfpts(:,2)-rnnew(n0,2)).^2+(surfpts(:,3)-rnnew(n0,3)).^2);
   newsurfNode = surfpts(point_id,:);

   % Split only the flag 0 node connected with the flag 6 node.
             % The velocity of flag 6 node is given to the new node to
             % compute new glide planes of new links in splitnode().
             if ~isempty(surfpts)
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
                    linksnew(end,6:8)=linksnew(connectivitynew(n0,2*ii),6:8);
                else %Post-force evaluation, rnnew nx6, [rn vel flag]
                    vntemp = rnnew(n1,4:6);
                    posvel = [newsurfNode, vntemp];
                    [rnnew,linksnew,connectivitynew,linksinconnectnew]=splitnode(rnnew,linksnew,connectivitynew,linksinconnectnew,n0,ii,posvel);
                    rnnew(end,end) = 6; %splitnode always gives flag 0, so we have to change it as 6 for surface node!
                    %testing beta
                    linksnew(end,6:8)=linksnew(connectivitynew(n0,2*ii),6:8);
                end
             end
end
