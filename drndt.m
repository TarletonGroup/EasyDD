function [vnvec,fn,fseg] = drndt(rnvec,flag,MU,NU,a,Ec,links,connectivity,...
    mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d)

% This needs to be an input/obtained from the surface nodes. This is a temporary fix for cuboid.
normals=[1 0 0;
          0 1 0;
          0 0 1];
tol=1e-6;

%unscramble rn
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

   % Virtual and fixed nodes have zero velocity
    if rn(p,4) == 7 || rn(p,4) == 67
        vn(p,:)=[0 0 0];
   % Surface nodes are confined to moving in the movement plane and surface plane.
    elseif rn(p,4) == 6
            connodes = [rn(links(links(:,1)==p,2),[1:3,end]); rn(links(links(:,2)==p,1),[1:3,end])];
            virtconnodes = connodes(connodes(:,4)==67,1:3);
            realconnodes = connodes(connodes(:,4)~=67,1:3);
            % If the surface node doesn't have a  single defined movement plane, make its velocity zero and continue.
            if size(realconnodes,1)>1 || isempty(realconnodes)
               vn(p,:)=[0 0 0];
               continue
            end
            slipplane=cross((rn(p,1:3)-realconnodes),vn(p,:));
            slipplane=slipplane/norm(slipplane);
            reps=ones(size(virtconnodes,1),3);
            reps=rn(p,1:3).*reps;
            vec=virtconnodes-reps;
            vec = sum(vec,1);
            vec=normalize(vec, 2);
            dotprods=normals*vec';
            surfplanes=normals(abs(dotprods)>tol,:);
            if isempty(surfplanes)
               rn(p,4)=0;
               continue
            end
            slipplane=slipplane.*ones(size(surfplanes));
            lines=cross(surfplanes,slipplane,2);
            lines = normalize(lines, 2);
            if size(lines,1)>3
                vn(p,:)=[0 0 0];
                continue
            elseif size(lines,1)==2
               if 1 - abs(lines(1,:) * lines(2,:)') < eps
                   lines=lines(1,:);
               else
                   vn(p,:)=[0 0 0];
                   continue
               end
            end
            vn(p,:)=(vn(p,:)*lines').*lines;
    end
end

%make up the last
vnvec=reshape(vn,length(vn(:,1))*3,1);
