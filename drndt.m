function [vnvec,fn,fseg] = drndt(rnvec,flag,MU,NU,a,Ec,links,connectivity,...
    mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d)

% unscramble rn
rn=[reshape(rnvec,length(rnvec)/3,3),flag];

%rn(:,1:3)

%nodal driving force
fseg=segforcevec(MU,NU,a,Ec,rn,links,0,vertices,...
    uhat,nc,xnodes,D,mx,mz,w,h,d);

%mobility function
[vn,fn]=feval(mobility,fseg,rn,links,connectivity,[],[]);

% fixed nodes (flag==7) are not allowed to move.
% flag == 6, are only allowed to move on the surface they live at.
% flag == 67, are virtual nodes that are not allowed to move.
for p=1:size(vn,1)

    if rn(p,4) == 7 || rn(p,4) == 67
        vn(p,:)=[0 0 0];
    end
    
    % Daniel Celis Garza 20200120: Surface node has zero velocity in the 
    % direction of the surface normal.
    if rn(p,4) == 6
        % Find the coordinates of the nodes connected to the surface node.
        connodes = [rn(links(links(:,1)==p,2),[1:3,end]);
                    rn(links(links(:,2)==p,1),[1:3,end])];
        % Out of the connected nodes find if any of them is external (67).
        virtual_connection = connodes(connodes(:,end)==67,1:3);
        % If the surface node is connected to an external node, we
        % eliminate the surface node's velocity component in that
        % direction.
        if ~isempty(virtual_connection)
            % If the node is connected to more than one external node we
            % find the mean line direction and use that vector to correct
            % the velocity.
            if size(virtual_connection,1)>1
                virtual_connection = mean(virtual_connection);
                virtual_connection = virtual_connection/norm(virtual_connection);
            end
            % Vector rejection eliminates the velocity components in the
            % direction of the virtual segment.
            vn(p,:) = vn(p,:) - dot(vn(p,:), virtual_connection)*virtual_connection;
        end
    end
    
end

%make up the last
vnvec=reshape(vn,length(vn(:,1))*3,1);
