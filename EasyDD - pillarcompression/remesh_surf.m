%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bruce Bromage
% Oxford Micromechanics Group
% Department of Materials, University of Oxford
% bruce.bromage@materials.ox.ac.uk
% December 2019
%
% Flags:
% 0 = free moving internal node
% 7 = fixed internal node (no burger's vector conservation)
% 6 = on free surface
% 65 = exited node that requires remeshing
% 67 = outside of domain, fixed.
%
% inhull function - checks whether node is inside or outside (convex!)
% domain, defined by vertices (surfaces are triangulated from vertices)
%
% NB, this is overloaded for both rn and rnnew arrays. When specifying
% flags, always use rn(:,end) or rnnew(:,end).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=...
    remesh_surf(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,vertices,P,fn)

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
   if rnnew(n,end)==65 && rnnew(n,3)<=0
       vec=[0,0,0];
       connumb=connectivitynew(n,1);
       for m=1:connumb
           vec=vec+rnnew(linksnew(connectivitynew(n,2*m),3-connectivitynew(n,2*m+1)),1:3)-rnnew(n,1:3);
       end
       vec=rnnew(n,3).*(vec/vec(1,3));
       rnnew(n,1:3)=rnnew(n,1:3)-vec;
       rnnew(n,end)=7;
       if any(isnan(rnnew(n,1:3)))
           fprintf('Error fixing node to back end. See remesh_surf line 62')
           pause;
       end
   end
end

%Create surface nodes for newly exited nodes
for n=1:L1
    n0=nodelist(n); % the node id which is currently considered.
    % Find outside nodes connected with inside node with flag 0 or 7
    if rnnew(n0,end)==0 || rnnew(n0,end)==7
        numNbrs=conlist(n,1); % the number of neighbor nodes

        for i=1:numNbrs
            ii=conlist(n,i+1); % connectionid for this connection
            linkid=connectivitynew(n0,2*ii); % the link id
            posinlink=connectivitynew(n0,2*ii+1);
            n1=linksnew(linkid,3-posinlink); % the neighbor node id
            if rnnew(n1,end)==65
                [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew] = gensurfnode2(Index,fn,P,rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,n0,n1,i,ii,vertices);
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

    %move exited node back onto surface
    [rnnew] = movetosurf(rnnew,linksnew,i,vertices);

    %estimate surface centroid nearest to point
    [~,index] = min( (P(:,1) - rnnew(i,1)).^2 + (P(:,2) - rnnew(i,2)).^2 + (P(:,3) - rnnew(i,3)).^2 );

    %find surface normal associated with centroid
    Index(i)=index;

    %move exited node to surface
    %[rnnew] = movetosurf(rnnew,linksnew,i,vertices);

    %extend far away
    rnnew = extend(rnnew,linksnew,i,index,fn);
    if any(any(isnan(rnnew)))
        disp('The new extended node is NaN! See Line 46 in remesh_surf.m');
        pause;
    end

end

%find connected surface nodes and merge them
i=1;
while i<size(rnnew,1)
    if rnnew(i,end)~=6
        i=i+1;
        continue
    end
    check=0;
    numcon=connectivitynew(i,1);
    for j=1:numcon
        connode=linksnew(connectivitynew(i,2*j),3-connectivitynew(i,2*j+1));
        if rnnew(connode,end)==6
            [rnnew,connectivitynew,linksnew,linksinconnectnew,fsegnew,~]=mergenodes(rnnew,connectivitynew,linksnew,linksinconnectnew,fsegnew,connode,i,1,1,1,1);
%             if nodeid~=0
%                 [rnnew] = movetosurf(rnnew,linksnew,nodeid,vertices);
%             end
            i=1;
            check=1;
            break
        end
    end
    if check==1
        continue
    end
    i=i+1;
end

%find surface nodes connected only with virtual nodes and repair them
%accordingly
L1=size(rnnew,1);
for i=1:L1 %length(rnnew)  ET updated- check!
    if rnnew(i,end)~=6
        continue;
    end
    test=0;
    for j=1:connectivitynew(i,1)
        %estimate surface centroid nearest to point
        if length(Index)<i || Index(i) == 0 || isnan(Index(i))  %update
            [~,index] = min( (P(:,1) - rnnew(i,1)).^2 + (P(:,2) - rnnew(i,2)).^2 + (P(:,3) - rnnew(i,3)).^2 );
            Index(i)=index;
        end
        if rnnew(linksnew(connectivitynew(i,2*j),(3-connectivitynew(i,2*j+1))),end)==0 || rnnew(linksnew(connectivitynew(i,2*j),(3-connectivitynew(i,2*j+1))),end)==7 %if the surface node is not connected to a virtual node : break M
            break;
%         elseif dot(fn(Index(i),:),rnnew(linksnew(connectivitynew(i,2*j),(3-connectivitynew(i,2*j+1))),1:3)-rnnew(i,1:3))<0
%             [rnnew,linksnew,connectivitynew,linksinconnectnew] = gensurfnode2(Index,fn,P,rnnew,linksnew,connectivitynew,linksinconnectnew,i,linksnew(connectivitynew(i,2*j),(3-connectivitynew(i,2*j+1))),j,j,vertices);
%             break
        else
            test=test+1;
        end
    end
    if test==connectivitynew(i,1) %all linked nodes are virtual (67 flag)
        %extend far away
        rnnew = extend(rnnew,linksnew,i,index,fn);
    end
end

i=1; %ensure surface nodes only connect to one virtual node
while i<size(rnnew,1)
    if rnnew(i,end)~=6
       i=i+1;
       continue
    end
    if connectivitynew(i,1)<=2
       i=i+1;
       continue
    end
    for j=1:connectivitynew(i,1)-1
        nid1=linksnew(connectivitynew(i,2*j),3-connectivitynew(i,2*j+1));
        if rnnew(nid1,end)~=67
           continue
        end
        for k=j+1:connectivitynew(i,1)
            nid2=linksnew(connectivitynew(i,2*k),3-connectivitynew(i,2*k+1));
            if rnnew(nid2,end)~=67
                continue
            end
            [rnnew,connectivitynew,linksnew,linksinconnectnew,fsegnew,~]=mergenodes(rnnew,connectivitynew,linksnew,linksinconnectnew,fsegnew,nid2,nid1,1,1,1,1);
            i=1;
            break
        end
        break
    end
    i=i+1;
end

for i=1:size(rnnew,1) %create surface nodes where necessary
    if rnnew(i,end)==67
       for j=1:connectivitynew(i,1)
          con_id=linksnew(connectivitynew(i,2*j),3-connectivitynew(i,2*j+1));
          if rnnew(con_id,end)==0 || rnnew(con_id,end)==7
              [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew] = gensurfnode2(Index,fn,P,rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,con_id,i,connectivitynew(i,1),j,vertices);
          end
       end
    end
end

end

function rnnew = extend(rnnew,~,rn_id,plane_id,fn)
    %% Extend along surface normal
    %N.B. This method should not be used with Fivel's displacement code

    p_inf=1e7;
    proj_vector=fn(plane_id,1:3);
    rnnew(rn_id,1:3) = rnnew(rn_id,1:3) + p_inf*proj_vector; %far away...
    rnnew(rn_id,end) = 67;   %flag node as virtual
end

function [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew] = gensurfnode2(~,~,~,rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,n0,n1,~,ii,vertices)
    faces = [1,3,4,2;                                             %Faces of cuboid as defined by vertices
             5,6,8,7;
             2,4,8,6;
             1,5,7,3;
             1,2,6,5;
             3,7,8,4];
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
                    fsegnew(end+1,:)=[0 0 0 0 0 0];
                else %Post-force evaluation, rnnew nx6, [rn vel flag]
                    vntemp = rnnew(n1,4:6);
                    posvel = [newsurfNode, vntemp];
                    [rnnew,linksnew,connectivitynew,linksinconnectnew]=splitnode(rnnew,linksnew,connectivitynew,linksinconnectnew,n0,ii,posvel);
                    rnnew(end,end) = 6; %splitnode always gives flag 0, so we have to change it as 6 for surface node!
                    %testing beta
                    linksnew(end,6:8)=linksnew(connectivitynew(n0,2*ii),6:8);
                    fsegnew(end+1,:)=[0 0 0 0 0 0];
                end
             end
end

function [rnnew] = movetosurf(rnnew,~,i,vertices)
    faces = [1,3,4,2;                                             %Faces of cuboid as defined by vertices
             5,6,8,7;
             2,4,8,6;
             1,5,7,3;
             1,2,6,5;
             3,7,8,4];
%     connodes = [rnnew(linksnew(linksnew(:,1)==i,2),[1:3,end]);rnnew(linksnew(linksnew(:,2)==i,1),[1:3,end])];
%     connodes = connodes(connodes(:,4)==67,1:3);
%     reps=ones(size(connodes,1),3);
%     reps=rnnew(i,1:3).*reps;
%     vec=connodes-reps;
%     vec = sum(vec,1);
%     vec=vec/norm(vec);
%     vec(abs(vec)~=max(abs(vec)))=0;  %convert mean direction into closeset surface normal (needs to be changed for non-cuboid)
%     vec(abs(vec)==max(abs(vec)))=1;
%
%     vec = [rnnew(linksnew(linksnew(:,1)==i,2),1:3)-rnnew(i,1:3);rnnew(linksnew(linksnew(:,2)==i,1),1:3)-rnnew(i,1:3)];
%     vec = sum(vec,1)/size(vec,1);
    vec = rnnew(i,4:6);
    if norm(vec)< eps
        fprintf('Error moving exited node back to surface. See movetosurf in remesh_surf')
        pause
    end
    vec=vec/norm(vec);
    line = [rnnew(i,1:3),vec];
    surfpts = intersectLineMesh3d(line, vertices, faces);
    [~,point_id] = min((surfpts(:,1)-rnnew(i,1)).^2+(surfpts(:,2)-rnnew(i,2)).^2+(surfpts(:,3)-rnnew(i,3)).^2);
    newsurfNode = surfpts(point_id,:);
    if isempty(newsurfNode)
        rnnew(i,end) = 0;
    else
        rnnew(i,1:3) = newsurfNode;
    end
%     rnnew(i,1:3) = rnnew(i,1:3)+vec;

end