function [rn,links,connectivity,linksinconnect,fseg,colliding_segments]=collision_basic(...
    rn,links,connectivity,linksinconnect,fseg,mindist,MU,NU,a,Ec,mobility,vertices,...
    uhat,nc,xnodes,D,mx,mz,w,h,d,floop,n1s1,n2s1,n1s2,n2s2,s1,s2,~,lmin)
%floop to know wich loop has to be run

colliding_segments=1;
mindist2=mindist*mindist;
lrn2=length(rn(1,:));
lrn3=lrn2-1; %lnr2-1 to count every column execpt the one with the flag
% tol is the error factor of the calculation
tol=1e-12;
remcon=0;
lmin2=lmin*lmin;

% if connectivity(n1s1,1)==2 && connectivity(n2s1,1)==2 && connectivity(n1s2,1)==2 && connectivity(n2s2,1)==2
%     alt1=connectivity(n1s1,[2 4]);
%     other1=alt1(alt1~=s1);
%     alt2=connectivity(n2s1,[2 4]);
%     other2=alt2(alt2~=s1);
%     alt3=connectivity(n1s2,[2 4]);
%     other3=alt3(alt3~=s2);
%     alt4=connectivity(n2s2,[2 4]);
%     other4=alt4(alt4~=s2);
%     if other1==other3 || other1==other4 || other2==other3 || other2==other4
%         remcon=1;
%     end
% end

% check for two links colliding
if remcon==1
    fprintf('Remeshing conflict detected. Cancelling collision\n')
elseif floop==1   %run loop1 only if floop=1, if floop=2 it means that only loop2 needs to be run.
    % First collision in computed with collisioncheckermexmarielle information : n1s1,n2s1,n1s2,n2s2,s1,s2
    
    [dist2,ddist2dt,L1,L2]=mindistcalcmex(rn(n1s1,1:lrn3),rn(n2s1,1:lrn3),rn(n1s2,1:lrn3),rn(n2s2,1:lrn3));
    collision_condition_is_met=((dist2<mindist2)&(ddist2dt<-tol))|(dist2<tol);
    % there are two conditions here the first condition handles non planar collisions
    % the second conditions looks for coplanar collisions
    if collision_condition_is_met
        % links are unconnected and colliding
        % identify the first node to be merged
        vec=rn(n2s1,1:3)-rn(n1s1,1:3);
        close_to_n1s1=((L1*L1*(vec*vec'))<lmin2);
        close_to_n2s1=(((1-L1)*(1-L1)*(vec*vec'))<lmin2);
        %                 small_link=((vec*vec')<(4*mindist2));
        %                 tiny_link=((vec*vec')<(mindist2));
        % if collision point is close to one of the existing nodes use that node
        if close_to_n1s1 && L1<=0.5 %&& connectivity(n1s1,1)<5 || connectivity(n2s1,1)<5 && small_link && ~tiny_link
            mergenode1=n1s1;
        elseif close_to_n2s1 %&& connectivity(n2s1,1)<5 || connectivity(n1s1,1)<5 && small_link && ~tiny_link
            mergenode1=n2s1;
            %                 elseif tiny_link
            %                     fprintf('Error detected in collision_basic. See Line 52\n')
            %                     pause
            % %                     colliding_segments=0;
        else
            spnode=n1s1;
            splitconnection=linksinconnect(s1,1); % linki=s1 M
            if close_to_n1s1
                L1=mindist/norm(vec);
            elseif close_to_n2s1
                L1=1-(mindist/norm(vec));
            end
            posvel=rn(n1s1,1:lrn3).*(1-L1)+rn(n2s1,1:lrn3).*L1;
            [rn,links,connectivity,linksinconnect]=splitnode(rn,links,connectivity,linksinconnect,spnode,splitconnection,posvel);
            mergenode1=length(rn(:,1)); %nodeid of mergenode1 M
            %linknew=length(links(:,1)); %linkid of linknew M
            %                     links(linknew,6:8)=links(s1,6:8); %glide plane M
            fseg=[fseg;zeros(1,6)];
        end
        
        % identify the second node to be merged
        vec=rn(n2s2,1:3)-rn(n1s2,1:3);
        close_to_n1s2=((L2*L2*(vec*vec'))<lmin2);
        close_to_n2s2=(((1-L2)*(1-L2)*(vec*vec'))<lmin2);
        %                 small_link=((vec*vec')<(4*mindist2));
        %                 tiny_link=((vec*vec')<(mindist2));
        % if collision point is close to one of the existing nodes use that node
        if close_to_n1s2 && L2<=0.5 %&& connectivity(n1s2,1)<5 || connectivity(n2s2,1)<5 && small_link && ~tiny_link
            mergenode2=n1s2;
        elseif close_to_n2s2 %&& connectivity(n2s2,1)<5 || connectivity(n1s2,1)<5 && small_link && ~tiny_link
            mergenode2=n2s2;
            %                 elseif tiny_link
            %                     fprintf('Error detected in collision_basic. See Line 83\n')
            %                     pause
            % %                     colliding_segments=0;
        else
            spnode=n1s2;
            splitconnection=linksinconnect(s2,1); %linkj=s2 M
            if close_to_n1s2
                L2=mindist/norm(vec);
            elseif close_to_n2s2
                L2=1-(mindist/norm(vec));
            end
            posvel=rn(n1s2,1:lrn3).*(1-L2)+rn(n2s2,1:lrn3).*L2;
            [rn,links,connectivity,linksinconnect]=splitnode(rn,links,connectivity,linksinconnect,spnode,splitconnection,posvel);
            mergenode2=length(rn(:,1));
            %linknew=length(links(:,1));
            %links(linknew,6:8)=links(s2,6:8);
            fseg=[fseg;zeros(1,6)];
        end
        
        if colliding_segments==1
            % merge the two colliding nodes
            %                 disp(sprintf('node %d and node %d are colliding by two line collision',mergenode1,mergenode2))
            collisionpoint=findcollisionpoint(mergenode1,mergenode2,rn,connectivity,links);
            rn(mergenode1,1:lrn2)=[collisionpoint 0 0 0 max(rn(mergenode1,lrn2),rn(mergenode2,lrn2)) ];
            [rn,connectivity,links,linksinconnect,fseg,mergednodeid]=mergenodes(rn,connectivity,links,linksinconnect,fseg,mergenode1,mergenode2,MU,NU,a,Ec);
            if mergednodeid>0
                for k=1:connectivity(mergednodeid,1)
                    linkid=connectivity(mergednodeid,2*k);
                    fseg(linkid,:)=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 lrn2]),links,linkid,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);
                    othernode=links(linkid,3-connectivity(mergednodeid,2*k+1)); % 3-connectivity(mergednodeid,2*k+1) = 1 or 2, it corresponds to the position of the other node of the link ( the one which is not mergenode ) M
                    clist=[connectivity(othernode,1) linspace(1,connectivity(othernode,1),connectivity(othernode,1))];
                    [rn(othernode,4:6),~]=feval(mobility,fseg,rn,links,connectivity,othernode,clist);
                end
                numbcon=connectivity(mergednodeid,1);
                conlist=[numbcon linspace(1,numbcon,numbcon)];
                [rn(mergednodeid,4:6),~]=feval(mobility,fseg,rn,links,connectivity,mergednodeid,conlist);
            end
        end
    end
else
    [dist2,ddist2dt,L1,~]=mindistcalcmex(rn(n1s1,1:lrn3),rn(n2s1,1:lrn3),rn(n1s2,1:lrn3),rn(n1s2,1:lrn3));
    collision_condition_is_met=(dist2<mindist2)&(ddist2dt<-tol);
    if collision_condition_is_met
        % identify the first node to be merged
        mergenode1=n1s2;
        % identify the second node to be merged
        vec=rn(n2s1,1:3)-rn(n1s1,1:3);
        close_to_n1s1=((L1*L1*(vec*vec'))<lmin2);
        close_to_n2s1=(((1-L1)*(1-L1)*(vec*vec'))<lmin2);
        if links(s2,1)==n1s1 || links(s2,2)==n1s1
            hingenode=n1s1;
        else
            hingenode=n2s1;
        end
        %                     small_link=((vec*vec')<(4*mindist2));
        %                     tiny_link=((vec*vec')<(mindist2));
        % if collision point is close to one of the existing nodes use that node
        if close_to_n1s1 && L1<=0.5%&& connectivity(n1s1,1)<5 || connectivity(n2s1,1)<5 && small_link && ~tiny_link
            if ~isequal(n1s1,hingenode)
                mergenode2=n1s1;
            else
                vec2=vec/norm(vec);
                vec2=lmin.*vec2;
                newpoint=rn(n1s1,1:3)+vec2;
                newdist=norm(rn(n2s1,1:3)-newpoint);
                if newdist<lmin
                    mergenode2=n2s1;
                else
                    L1=1-(newdist/norm(vec));
                    spnode=n1s1;
                    splitconnection=linksinconnect(s1,1);
                    posvel=rn(n1s1,1:lrn3).*(1-L1)+rn(n2s1,1:lrn3).*L1;
                    [rn,links,connectivity,linksinconnect]=splitnode(rn,links,connectivity,linksinconnect,spnode,splitconnection,posvel);
                    mergenode2=length(rn(:,1));
                    newlink=length(links(:,1));
                    links(newlink,6:8)=links(s1,6:8);
                    fseg=[fseg;zeros(1,6)];
                end
            end
        elseif close_to_n2s1 %&& connectivity(n2s1,1)<5 || connectivity(n1s1,1)<5 && small_link && ~tiny_link
            if ~isequal(n2s1,hingenode)
                mergenode2=n2s1;
            else
                vec2=vec/norm(vec);
                vec2=lmin.*vec2;
                newpoint=rn(n2s1,1:3)-vec2;
                newdist=norm(rn(n1s1,1:3)-newpoint);
                if newdist<lmin
                   mergenode2=n1s1;
                else
                    L1=newdist/norm(vec);
                    spnode=n1s1;
                    splitconnection=linksinconnect(s1,1);
                    posvel=rn(n1s1,1:lrn3).*(1-L1)+rn(n2s1,1:lrn3).*L1;
                    [rn,links,connectivity,linksinconnect]=splitnode(rn,links,connectivity,linksinconnect,spnode,splitconnection,posvel);
                    mergenode2=length(rn(:,1));
                    newlink=length(links(:,1));
                    links(newlink,6:8)=links(s1,6:8);
                    fseg=[fseg;zeros(1,6)]; 
                end
            end
            %                     elseif tiny_link
            %                         fprintf('Error detected in collision_basic. See Line 142\n')
            %                         pause
        else
            spnode=n1s1;
            splitconnection=linksinconnect(s1,1);
            posvel=rn(n1s1,1:lrn3).*(1-L1)+rn(n2s1,1:lrn3).*L1;
            [rn,links,connectivity,linksinconnect]=splitnode(rn,links,connectivity,linksinconnect,spnode,splitconnection,posvel);
            mergenode2=length(rn(:,1));
            newlink=length(links(:,1));
            links(newlink,6:8)=links(s1,6:8);
            fseg=[fseg;zeros(1,6)];
        end
        %merge the two nodes
        %                     disp(sprintf('node %d and node %d are colliding by hinge condition',mergenode2,mergenode1))
        collisionpoint=findcollisionpoint(mergenode1,mergenode2,rn,connectivity,links);
        rn(mergenode1,1:lrn2)=[collisionpoint 0 0 0 max(rn(mergenode1,lrn2),rn(mergenode2,lrn2)) ];
        [rn,connectivity,links,linksinconnect,fseg,mergednodeid]=mergenodes(rn,connectivity,links,linksinconnect,fseg,mergenode1,mergenode2,MU,NU,a,Ec);
        if mergednodeid>0
            for k=1:connectivity(mergednodeid,1)
                linkid=connectivity(mergednodeid,2*k);
                fseg(linkid,:)=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 lrn2]),links,linkid,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);
                othernode=links(linkid,3-connectivity(mergednodeid,2*k+1));
                clist=[connectivity(othernode,1) linspace(1,connectivity(othernode,1),connectivity(othernode,1))];
                [rn(othernode,4:6),~]=feval(mobility,fseg,rn,links,connectivity,othernode,clist);
            end
            numbcon=connectivity(mergednodeid,1);
            conlist=[numbcon linspace(1,numbcon,numbcon)];
            [rn(mergednodeid,4:6),~]=feval(mobility,fseg,rn,links,connectivity,mergednodeid,conlist);
        end
    end
end
function collisionpoint=findcollisionpoint(mergenode1,mergenode2,rn,connectivity,links)
% this subroutine finds the collision point of two nodes given that there are strict glide plane constraints
eps=1e-12;
newplanecondition=0.875;
p1=rn(mergenode1,1:3);
p2=rn(mergenode2,1:3);
Nmat=zeros(3,3);
Nsize=0;
vector=zeros(3,1);
s=size(rn,2);
if rn(mergenode1,s)==7
    collisionpoint=rn(mergenode1,1:3);
    return;
elseif rn(mergenode2,s)==7
    collisionpoint=rn(mergenode2,1:3);
    return;
end

for i=1:connectivity(mergenode1,1)
    if Nsize<3
        linkid=connectivity(mergenode1,2*i);
        connode=links(connectivity(mergenode1,2*i),3-connectivity(mergenode1,2*i+1));
        rt=rn(mergenode1,1:3)-rn(connode,1:3);
        L=norm(rt);
        linedir=rt./L;
        n1=cross(linedir,links(linkid,3:5));
        n2=cross(linedir,rn(mergenode1,4:6));
        
        if n1*n1'>eps
            plane=n1./norm(n1);
        elseif n2*n2'>eps
            plane=n2./norm(n2);
        end
        
        if ((n1*n1'>eps)||(n2*n2'>eps))
            if Nsize==0
                conditionismet = 1;
            elseif Nsize==1
                conditionismet = ((Nmat(1,:)*plane')^2 < newplanecondition*newplanecondition);
            else
                detN=det([Nmat(1:2,:);plane]);
                conditionismet = detN*detN > (1-newplanecondition)^4;
            end
            if conditionismet
                Nsize=Nsize+1;
                Nmat(Nsize,:)=plane;
                vector(Nsize)=plane*p1';
            end
        end
    end
end

for i=1:connectivity(mergenode2,1)
    if Nsize<3
        linkid=connectivity(mergenode2,2*i);
        connode=links(connectivity(mergenode2,2*i),3-connectivity(mergenode2,2*i+1));
        rt=rn(mergenode2,1:3)-rn(connode,1:3);
        L=norm(rt);
        linedir=rt./L;
        n1=cross(linedir,links(linkid,3:5));
        n2=cross(linedir,rn(mergenode2,4:6));
        
        if n1*n1'>eps
            plane=n1./norm(n1);
        elseif n2*n2'>eps
            plane=n2./norm(n2);
        end
        if ((n1*n1'>eps)||(n2*n2'>eps))
            if Nsize==1
                conditionismet = ((Nmat(1,:)*plane')^2 < newplanecondition*newplanecondition);
            else
                detN=det([Nmat(1:2,:);plane]);
                conditionismet = detN*detN > (1-newplanecondition)^4;
            end
            
            if conditionismet
                Nsize=Nsize+1;
                Nmat(Nsize,:)=plane;
                vector(Nsize)=plane*p2';
            end
        end
    end
end

Matrix=[eye(3) Nmat(1:Nsize,:)';Nmat(1:Nsize,:) zeros(Nsize,Nsize)];
V=[(rn(mergenode1,1:3)'+rn(mergenode2,1:3)')./2; vector(1:Nsize)];
res=Matrix\V;
collisionpoint=res(1:3)';
