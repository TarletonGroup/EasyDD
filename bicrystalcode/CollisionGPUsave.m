%Marielle Thibault
% this subroutine goes through the existing links and checks for collisions
% it first checks through unconnected links
% it then checks for hinges that are coming within the mininum distance
% Same role as collision but is use GPU
% The main steps in CollisionGPU.m the new collision function are:
% -	Creating a list with all possibilities = inputmin
% -	Running mindistcalcGPU1.m : algebraic computations
% -	Running mindistcalcGPU2.m : it gives for each possibility the same results as mindistcalc.m that is to say L1,L2,dist2,ddist2dt. 
% -	Check for every possibility if condition_collision_is_met and if it is the case do the remeshing. 

% Fonctions created for collisionGPU.m :
% - CreateInputMex : creation of inputmin
% - MindistcalcGPU1.m
% - MindistcalcGPU2.m
% - CollisionGPU.m works only if it is CollisionCheckerMexMarielle.c


function [rn,links,connectivity,linksinconnect,fseg]=collisionGPUsave(rn,links,connectivity,linksinconnect,fseg,mindist,MU,NU,a,Ec,mobility,vertices,...
    uhat,nc,xnodes,D,mx,mz,w,h,d,floop,n1s1,n2s1,n1s2,n2s2,s1,s2,segpair)
%floop to know wich loop has to be run 

mindist2=mindist*mindist;
lrn2=length(rn(1,:));
lrn3=lrn2-1; %lnr2-1 to count every column execpt the one with the flag
% eps is the error factor of the calculation
eps=1e-12;

% check for two links colliding

if floop==1   %run loop1 only if floop=1, if floop=2 it means that only loop2 needs to be run. 
    
%First collision in computed with collisioncheckermexmarielle information : n1s1,n2s1,n1s2,n2s2,s1,s2

[dist2,ddist2dt,L1,L2]=mindistcalcmex(rn(n1s1,1:lrn3),rn(n2s1,1:lrn3),rn(n1s2,1:lrn3),rn(n2s2,1:lrn3));
            collision_condition_is_met=((dist2<mindist2)&(ddist2dt<-eps))|(dist2<eps);
            % there are two conditions here the first condition handles non planar collisions
            % the second conditions looks for coplanar collisions
            if collision_condition_is_met
                % links are unconnected and colliding
                % identify the first node to be merged
                vec=rn(n1s1,1:3)-rn(n2s1,1:3);
                close_to_n1s1=((L1*L1*(vec*vec'))<mindist2);
                close_to_n2s1=(((1-L1)*(1-L1)*(vec*vec'))<mindist2);
                % if collision point is close to one of the existing nodes use that node
                if close_to_n1s1
                    mergenode1=n1s1;
                elseif close_to_n2s1
                    mergenode1=n2s1;
                else
                    spnode=n1s1;
                    splitconnection=linksinconnect(s1,1); % linki=s1 M
                    posvel=rn(n1s1,1:lrn3).*(1-L1)+rn(n2s1,1:lrn3).*L1;
                    [rn,links,connectivity,linksinconnect]=splitnode(rn,links,connectivity,linksinconnect,spnode,splitconnection,posvel);
                    mergenode1=length(rn(:,1)); %nodeid of mergenode1 M
                    linknew=length(links(:,1)); %linkid of linknew M
                    links(linknew,6:8)=links(s1,6:8); %glide plane M 
                    fseg=[fseg;zeros(1,6)];
                end

                % identify the second node to be merged
                vec=rn(n1s2,1:3)-rn(n2s2,1:3);
                close_to_n1s2=((L2*L2*(vec*vec'))<mindist2);
                close_to_n2s2=(((1-L2)*(1-L2)*(vec*vec'))<mindist2);
                % if collision point is close to one of the existing nodes use that node
                if close_to_n1s2 
                    mergenode2=n1s2;
                elseif close_to_n2s2
                    mergenode2=n2s2;
                else
                    spnode=n1s2;
                    splitconnection=linksinconnect(s2,1); %linkj=s2 M
                    posvel=rn(n1s2,1:lrn3).*(1-L2)+rn(n2s2,1:lrn3).*L2;
                    [rn,links,connectivity,linksinconnect]=splitnode(rn,links,connectivity,linksinconnect,spnode,splitconnection,posvel);
                    mergenode2=length(rn(:,1));
                    linknew=length(links(:,1));
                    links(linknew,6:8)=links(s2,6:8);
                    fseg=[fseg;zeros(1,6)];
                end
                % merge the two colliding nodes
                disp(sprintf('node %d and node %d are colliding by two line collision',mergenode1,mergenode2))
                collisionpoint=findcollisionpoint(mergenode1,mergenode2,rn,connectivity,links);
                rn(mergenode1,1:lrn2)=[collisionpoint 0 0 0 max(rn(mergenode1,lrn2),rn(mergenode2,lrn2)) ];
                [rn,connectivity,links,linksinconnect,fseg,mergednodeid]=mergenodes(rn,connectivity,links,linksinconnect,fseg,mergenode1,mergenode2,MU,NU,a,Ec);
                if mergednodeid>0
                    for k=1:connectivity(mergednodeid,1)
                        linkid=connectivity(mergednodeid,2*k);
                        fseg(linkid,:)=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 lrn2]),links,linkid,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);
                        othernode=links(linkid,3-connectivity(mergednodeid,2*k+1)); % 3-connectivity(mergednodeid,2*k+1) = 1 or 2, it corresponds to the position of the other node of the link ( the one which is not mergenode ) M 
                        clist=[connectivity(othernode,1) linspace(1,connectivity(othernode,1),connectivity(othernode,1))];
                        [rn(othernode,4:6),fntmp]=feval(mobility,fseg,rn,links,connectivity,othernode,clist);
                    end
                    numbcon=connectivity(mergednodeid,1);
                    conlist=[numbcon linspace(1,numbcon,numbcon)];
                    [rn(mergednodeid,4:6),fntmp]=feval(mobility,fseg,rn,links,connectivity,mergednodeid,conlist);
                end
            end  
           
%For the other collision now... 
    
%pre-heating of the computer...
A=2;
A=gpuArray(A);
A=gather(A);
           
    
            %inputmin with MEX file
            [r11x,r11y,r11z,r21x,r21y,r21z,r12x,r12y,r12z,r22x,r22y,r22z,v11x,v11y,v11z,v21x,v21y,v21z,v12x,v12y,v12z,v22x,v22y,v22z,...
            in1s1,in2s1,in1s2,in2s2,ii,ij,p]=CreateInputMex(rn(:,1),rn(:,2),rn(:,3),rn(:,4),rn(:,5),rn(:,6),rn(:,end),links(:,1),links(:,2));
        
            inputmin=[r11x(1:p),r11y(1:p),r11z(1:p),r21x(1:p),r21y(1:p),r21z(1:p),r12x(1:p),r12y(1:p),r12z(1:p),r22x(1:p),r22y(1:p),r22z(1:p),...
            v11x(1:p),v11y(1:p),v11z(1:p),v21x(1:p),v21y(1:p),v21z(1:p),v12x(1:p),v12y(1:p),v12z(1:p),v22x(1:p),v22y(1:p),v22z(1:p),...
            in1s1(1:p),in2s1(1:p),in1s2(1:p),in2s2(1:p),ii(1:p),ij(1:p)];        
          
            %segpair : subscript of the segment pair that has been collided with collisioncheckemexmarielle data
            R=zeros(p-segpair,13,'gpuArray'); %GPU array mindistcalcpart1 
            T=zeros(p-segpair,4,'gpuArray'); %GPU array mindistcalcpart2
            inputmin=gpuArray(inputmin); %GPU array input for mindistcalcpart1
            segpair=segpair+1;
             
            %mindistpart1 for all the segment pairs remaining from :segpair+1 to p
            [R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6),R(:,7),R(:,8),R(:,9),R(:,10),R(:,11),R(:,12),R(:,13)]...
            =arrayfun(@mindistcalcGPU1,inputmin(segpair:p,1),inputmin(segpair:p,2),inputmin(segpair:p,3),inputmin(segpair:p,4),inputmin(segpair:p,5),inputmin(segpair:p,6),inputmin(segpair:p,7),inputmin(segpair:p,8),inputmin(segpair:p,9),inputmin(segpair:p,10),inputmin(segpair:p,11),inputmin(segpair:p,12));

            %mindistpart2 for all the segment pairs remaining 
            [T(:,1),T(:,2),T(:,3),T(:,4)] =arrayfun(@mindistcalcGPU2,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6),R(:,7),R(:,8),R(:,9),R(:,10),R(:,11),R(:,12),R(:,13),...
            inputmin(segpair:p,1),inputmin(segpair:p,2),inputmin(segpair:p,3),inputmin(segpair:p,4),inputmin(segpair:p,5),inputmin(segpair:p,6),inputmin(segpair:p,7),inputmin(segpair:p,8),inputmin(segpair:p,9),inputmin(segpair:p,10),inputmin(segpair:p,11),inputmin(segpair:p,12),...
            inputmin(segpair:p,13),inputmin(segpair:p,14),inputmin(segpair:p,15),inputmin(segpair:p,16),inputmin(segpair:p,17),inputmin(segpair:p,18),inputmin(segpair:p,19),inputmin(segpair:p,20),inputmin(segpair:p,21),inputmin(segpair:p,22),inputmin(segpair:p,23),inputmin(segpair:p,24));
            %END
            T=gather(T);
            inputmin=gather(inputmin(:,25:30));
           
            kk=1;
            stop=segpair-1;
            while kk<=(p-stop)
                
               dist2=T(kk,3); %M
               ddist2dt=T(kk,4); %M
               collision_condition_is_met=((dist2<mindist2)&(ddist2dt<-eps))|(dist2<eps);
                % there are two conditions here the first condition handles non planar collisions
                % the second conditions looks for coplanar collisions
                if collision_condition_is_met
                    
                    L1=T(kk,1); %M
                    L2=T(kk,2); %M
                    n1s1=inputmin(kk+stop,1); %M
                    n2s1=inputmin(kk+stop,2); %M
                    n1s2=inputmin(kk+stop,3); %M
                    n2s2=inputmin(kk+stop,4); %M
                    s1=inputmin(kk+stop,5); %M
                    s2=inputmin(kk+stop,6); %M
                    
                    % links are unconnected and colliding
                    % identify the first node to be merged
                    vec=rn(n1s1,1:3)-rn(n2s1,1:3);
                    close_to_n1s1=((L1*L1*(vec*vec'))<mindist2);
                    close_to_n2s1=(((1-L1)*(1-L1)*(vec*vec'))<mindist2);
                    % if collision point is close to one of the existing nodes use that node
                    if close_to_n1s1
                        mergenode1=n1s1;
                    elseif close_to_n2s1
                        mergenode1=n2s1;
                    else
                        spnode=n1s1;
                        splitconnection=linksinconnect(s1,1); % linki=s1 M
                        posvel=rn(n1s1,1:6).*(1-L1)+rn(n2s1,1:6).*L1;
                        [rn,links,connectivity,linksinconnect]=splitnode(rn,links,connectivity,linksinconnect,spnode,splitconnection,posvel);
                        mergenode1=length(rn(:,1)); %nodeid of mergenode1 M
                        linknew=length(links(:,1)); %linkid of linknew M
                        links(linknew,6:8)=links(s1,6:8); %glide plane M 
                        fseg=[fseg;zeros(1,6)];
                    end
                    % identify the second node to be merged
                    vec=rn(n1s2,1:3)-rn(n2s2,1:3);
                    close_to_n1s2=((L2*L2*(vec*vec'))<mindist2);
                    close_to_n2s2=(((1-L2)*(1-L2)*(vec*vec'))<mindist2);
                    % if collision point is close to one of the existing nodes use that node
                    
                    if close_to_n1s2 
                        mergenode2=n1s2;
                    elseif close_to_n2s2
                        mergenode2=n2s2;
                    else
                        spnode=n1s2;         
                        splitconnection=linksinconnect(s2,1); %linkj=s2 M BUGG LINKSINCONNECT
                        posvel=rn(n1s2,1:6).*(1-L2)+rn(n2s2,1:6).*L2;
                        [rn,links,connectivity,linksinconnect]=splitnode(rn,links,connectivity,linksinconnect,spnode,splitconnection,posvel);
                        mergenode2=length(rn(:,1));
                        linknew=length(links(:,1));
                        links(linknew,6:8)=links(s2,6:8);
                        fseg=[fseg;zeros(1,6)];
                    end
                    
                     % merge the two colliding nodes
                    disp(sprintf('node %d and node %d are colliding by two line collision',mergenode1,mergenode2))
                    collisionpoint=findcollisionpoint(mergenode1,mergenode2,rn,connectivity,links);
                    rn(mergenode1,1:lrn2)=[collisionpoint 0 0 0 max(rn(mergenode1,lrn2),rn(mergenode2,lrn2)) ];
                    [rn,connectivity,links,linksinconnect,fseg,mergednodeid]=mergenodes(rn,connectivity,links,linksinconnect,fseg,mergenode1,mergenode2,MU,NU,a,Ec);
                    if mergednodeid>0
                        for k=1:connectivity(mergednodeid,1)
                            linkid=connectivity(mergednodeid,2*k);
                            fseg(linkid,:)=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 lrn2]),links,linkid,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);
                            othernode=links(linkid,3-connectivity(mergednodeid,2*k+1)); % 3-connectivity(mergednodeid,2*k+1) = 1 or 2, it corresponds to the position of the other node of the link ( the one which is not mergenode ) M 
                            clist=[connectivity(othernode,1) linspace(1,connectivity(othernode,1),connectivity(othernode,1))];
                            [rn(othernode,4:6),fntmp]=feval(mobility,fseg,rn,links,connectivity,othernode,clist);
                        end
                        numbcon=connectivity(mergednodeid,1);
                        conlist=[numbcon linspace(1,numbcon,numbcon)];
                        [rn(mergednodeid,4:6),fntmp]=feval(mobility,fseg,rn,links,connectivity,mergednodeid,conlist);
                    end
                    
                    %inputmin with MEX file
                    [r11x,r11y,r11z,r21x,r21y,r21z,r12x,r12y,r12z,r22x,r22y,r22z,v11x,v11y,v11z,v21x,v21y,v21z,v12x,v12y,v12z,v22x,v22y,v22z,...
                    in1s1,in2s1,in1s2,in2s2,ii,ij,p]=CreateInputMex(rn(:,1),rn(:,2),rn(:,3),rn(:,4),rn(:,5),rn(:,6),rn(:,end),links(:,1),links(:,2));

                    inputmin=[r11x(1:p),r11y(1:p),r11z(1:p),r21x(1:p),r21y(1:p),r21z(1:p),r12x(1:p),r12y(1:p),r12z(1:p),r22x(1:p),r22y(1:p),r22z(1:p),...
                    v11x(1:p),v11y(1:p),v11z(1:p),v21x(1:p),v21y(1:p),v21z(1:p),v12x(1:p),v12y(1:p),v12z(1:p),v22x(1:p),v22y(1:p),v22z(1:p),...
                    in1s1(1:p),in2s1(1:p),in1s2(1:p),in2s2(1:p),ii(1:p),ij(1:p)];        

                    R=zeros(p-kk,13,'gpuArray'); %GPU array mindistcalcpart1
                    T=zeros(p-kk,4,'gpuArray'); %GPU array mindistcalcpart2
                    inputmin=gpuArray(inputmin); %GPU array input for mindistcalcpart1 
                    
                    kk=kk+1;
                    %mindistpart1
                    [R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6),R(:,7),R(:,8),R(:,9),R(:,10),R(:,11),R(:,12),R(:,13)]...
                    =arrayfun(@mindistcalcGPU1,inputmin(kk:p,1),inputmin(kk:p,2),inputmin(kk:p,3),inputmin(kk:p,4),inputmin(kk:p,5),inputmin(kk:p,6),inputmin(kk:p,7),inputmin(kk:p,8),inputmin(kk:p,9),inputmin(kk:p,10),inputmin(kk:p,11),inputmin(kk:p,12));
                    %mindistpart2
                    [T(:,1),T(:,2),T(:,3),T(:,4)] =arrayfun(@mindistcalcGPU2,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6),R(:,7),R(:,8),R(:,9),R(:,10),R(:,11),R(:,12),R(:,13),...
                    inputmin(kk:p,1),inputmin(kk:p,2),inputmin(kk:p,3),inputmin(kk:p,4),inputmin(kk:p,5),inputmin(kk:p,6),inputmin(kk:p,7),inputmin(kk:p,8),inputmin(kk:p,9),inputmin(kk:p,10),inputmin(kk:p,11),inputmin(kk:p,12),...
                    inputmin(kk:p,13),inputmin(kk:p,14),inputmin(kk:p,15),inputmin(kk:p,16),inputmin(kk:p,17),inputmin(kk:p,18),inputmin(kk:p,19),inputmin(kk:p,20),inputmin(kk:p,21),inputmin(kk:p,22),inputmin(kk:p,23),inputmin(kk:p,24));
                    %end
                    
                    T=gather(T);
                    inputmin=gather(inputmin(:,25:30));
                    stop=kk-1; %for subscript in inputmin M 
                    kk=0;

                end
            kk=kk+1;   
            end

end        
        
% check for a hinge condition
i=1;
while i<=length(rn(:,1))
    %ignore virtual nodes - FF
    if rn(i,end)==67
        i=i+1;
        continue;
    end
    j=1;
    while j<=connectivity(i,1)
        nodenoti=links(connectivity(i,2*j),3-connectivity(i,2*j+1));
        k=1;
        while k<=connectivity(i,1)
            linkid=connectivity(i,2*k);
            % if node is on the link do not check for collision
            if j~=k
                % identify the nodes on the link
                n1s1=links(linkid,1);
                n2s1=links(linkid,2);
                %[dist2,ddist2dt,L1,L2]=mindistcalc(rn(n1s1,1:lrn3),rn(n2s1,1:lrn3),rn(nodenoti,1:lrn3),rn(nodenoti,1:lrn3));
                [dist2,ddist2dt,L1,L2]=mindistcalcmex(rn(n1s1,1:lrn3),rn(n2s1,1:lrn3),rn(nodenoti,1:lrn3),rn(nodenoti,1:lrn3));
                %dist2
                %ddist2dt
                collision_condition_is_met=(dist2<mindist2)&(ddist2dt<-eps);
                if collision_condition_is_met
                    % identify the first node to be merged
                    mergenode1=nodenoti;
                    % identify the second node to be merged
                    vec=rn(n1s1,1:3)-rn(n2s1,1:3);
                    close_to_n1s1=((L1*L1*(vec*vec'))<mindist2);
                    close_to_n2s1=(((1-L1)*(1-L1)*(vec*vec'))<mindist2);
                    % if collision point is close to one of the existing nodes use that node
                    if close_to_n1s1
                        mergenode2=n1s1;
                    elseif close_to_n2s1
                        mergenode2=n2s1;
                    else
                        spnode=n1s1;
                        splitconnection=linksinconnect(linkid,1);
                        posvel=rn(n1s1,1:lrn3).*(1-L1)+rn(n2s1,1:lrn3).*L1;
                        [rn,links,connectivity,linksinconnect]=splitnode(rn,links,connectivity,linksinconnect,spnode,splitconnection,posvel);
                        mergenode2=length(rn(:,1));
                        newlink=length(links(:,1));
                        links(newlink,6:8)=links(linkid,6:8);
                        fseg=[fseg;zeros(1,6)];
                    end
                    %merge the two nodes
                    disp(sprintf('node %d and node %d are colliding by hinge condition',mergenode2,mergenode1))
                    collisionpoint=findcollisionpoint(mergenode1,mergenode2,rn,connectivity,links);
                    rn(mergenode1,1:lrn2)=[collisionpoint 0 0 0 max(rn(mergenode1,lrn2),rn(mergenode2,lrn2)) ];
                    [rn,connectivity,links,linksinconnect,fseg,mergednodeid]=mergenodes(rn,connectivity,links,linksinconnect,fseg,mergenode1,mergenode2,MU,NU,a,Ec);
                    if mergednodeid>0
                        for k=1:connectivity(mergednodeid,1)
                            linkid=connectivity(mergednodeid,2*k);
                            fseg(linkid,:)=segforcevec(MU,NU,a,Ec,rn(:,[1 2 3 lrn2]),links,linkid,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);
                            othernode=links(linkid,3-connectivity(mergednodeid,2*k+1));
                            clist=[connectivity(othernode,1) linspace(1,connectivity(othernode,1),connectivity(othernode,1))];
                            [rn(othernode,4:6),fntmp]=feval(mobility,fseg,rn,links,connectivity,othernode,clist);
                        end
                        numbcon=connectivity(mergednodeid,1);
                        conlist=[numbcon linspace(1,numbcon,numbcon)];
                        [rn(mergednodeid,4:6),fntmp]=feval(mobility,fseg,rn,links,connectivity,mergednodeid,conlist);
                    end 
                    %there has been a connectivity change in node i start the search through node i's connections from the beginning
                    if i>size(rn,1)
                        % this is a rare but possible case.
                        % for this condition to be satisfied the last node was being checked for closed hinge condition and it merged with another node
                        % since this was the last node being checked exit the function
                        return;
                    else
                        j=0;
                        k=connectivity(i,1);
                    end
                end
            end
            k=k+1;
        end
        j=j+1;
    end
    i=i+1;
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
         
         if ((n1*n1'>eps)|(n2*n2'>eps))
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