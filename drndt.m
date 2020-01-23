function [vnvec,fn,fseg] = drndt(rnvec,flag,MU,NU,a,Ec,links,connectivity,...
    mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d)

%unscramble rn
rn=[reshape(rnvec,length(rnvec)/3,3),flag];

%rn(:,1:3)

%nodal driving force
fseg=segforcevec(MU,NU,a,Ec,rn,links,0,vertices,...
    uhat,nc,xnodes,D,mx,mz,w,h,d);

%mobility function
[vn,fn]=feval(mobility,fseg,rn,links,connectivity,[],[]);

%fixed nodes (flag==7) are not allowed to move.
% flag == 6, are only allowed to move on the surface they live at.
for p=1:size(vn,1)

    if rn(p,4) == 7 || rn(p,4) == 67

        vn(p,:)=[0 0 0];

    end
    
    % Daniel Celis Garza 20200120: Surface node has zero velocity in the 
    % direction of the surface normal. Works for a cuboid.
    if rn(p,4) == 6
        connodes = [rn(links(links(:,1)==p,2),[1:3,end]);
                    rn(links(links(:,2)==p,1),[1:3,end])];
        % we want to find the external node/nodes (67) connected to the surface
        % node because it will tell us which surface it should be on.
        % We only do this if the surface node is connected to a virtual
        % node because there are at times surface nodes that are not on
        % the surface as a result of remeshing. We want those to keep
        % moving.
        virtual_connection = connodes(connodes(:,end)==67,1:3);
        if ~isempty(virtual_connection)
            if size(virtual_connection,1)>1
                virtual_connection = mean(virtual_connection);
            end
            [~, idx] = max(abs(virtual_connection));
            vn(p,idx) = 0;
        end
    end
    
end

%make up the last
vnvec=reshape(vn,length(vn(:,1))*3,1);

%vn
%disp(sprintf('t=%20.10e vn(1,:)=(%20.10e,%20.10e,%20.10e)',t,vn(1,1),vn(1,2),vn(1,3)));
%pause(0.2)
%if(t>0)
%   pause
%end