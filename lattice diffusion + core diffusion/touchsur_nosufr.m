function [rn,links,connectivity,linksinconnect,fseg]=touchsur(rn,links,connectivity,linksinconnect,fseg,mindist,MU,NU,a,Ec,mobility,appliedstress,surpos,lable);
    lrn2=length(rn(1,:));
    lrn3=lrn2-1;
    %llink=size(links,1);
    i=1;
    while i<=size(links,1)
        n1=links(i,1); %%%bug link increase so quickly, less than i of last step change to while.
        n2=links(i,2);
        %---------------------------------------------
        %  free surface in X Y Z direction
        %----------------------------------------------
        if abs(rn(n1,lable))>surpos |  abs(rn(n2,lable))>surpos
          if(rn(n2,lable)~=rn(n1,lable))
            if (rn(n1,lable)>0) | (rn(n2,lable)>0)
                L1=( surpos-rn(n1,lable) ) / ( rn(n2,lable)-rn(n1,lable) );
            elseif (rn(n1,lable)<0) | (rn(n2,lable)<0)
                L1=( -surpos-rn(n1,lable) ) / ( rn(n2,lable)-rn(n1,lable) );
            end
            posvel=rn(n1,1:lrn3).*(1-L1)+rn(n2,1:lrn3).*L1;
          end
        end
        
        if (abs(rn(n1,lable))>surpos) & (abs( rn(n2,lable))<=surpos )
            if connectivity(n1,1)>1
                spnode=n1;
                splitconnection=linksinconnect(i,1);
                
                [rn,links,connectivity,linksinconnect]=separatenode(rn,links,connectivity,linksinconnect,spnode,splitconnection,posvel); 
                if ( posvel(lable)==rn(n2,lable) )
                    mergenode1=n2;
                    mergenode2=length(rn(:,1));
                    [rn,connectivity,links,linksinconnect,fseg,mergednodeid]=mergenodes(rn,connectivity,links,linksinconnect,fseg,mergenode1,mergenode2,MU,NU,a,Ec);    
                    %[links,connectivity,linksinconnect,fseg]=removelink(links,connectivity,linksinconnect,fseg,i);
                end
                
            else
                rn(n1,1:lrn3)=posvel;
                if ( posvel(lable)==rn(n2,lable) )
                    mergenode1=n2;
                    mergenode2=n1;
                    [rn,connectivity,links,linksinconnect,fseg,mergednodeid]=mergenodes(rn,connectivity,links,linksinconnect,fseg,mergenode1,mergenode2,MU,NU,a,Ec);    
                    %[links,connectivity,linksinconnect,fseg]=removelink(links,connectivity,linksinconnect,fseg,i);    
                end
            end
        elseif  (abs(rn(n2,lable))>surpos) & (abs( rn(n1,lable))<=surpos )
             if connectivity(n2,1)>1
                spnode=n2;
                splitconnection=linksinconnect(i,2);
                [rn,links,connectivity,linksinconnect]=separatenode(rn,links,connectivity,linksinconnect,spnode,splitconnection,posvel); 
                if ( posvel(lable)==rn(n1,lable) )
                    mergenode1=n1;
                    mergenode2=length(rn(:,1));
                    [rn,connectivity,links,linksinconnect,fseg,mergednodeid]=mergenodes(rn,connectivity,links,linksinconnect,fseg,mergenode1,mergenode2,MU,NU,a,Ec);    
                    %[links,connectivity,linksinconnect,fseg]=removelink(links,connectivity,linksinconnect,fseg,i);    
                end
            else
                rn(n2,1:lrn3)=posvel;
                if ( posvel(lable)==rn(n1,lable) )
                    mergenode1=n1;
                    mergenode2=n2;
                    [rn,connectivity,links,linksinconnect,fseg,mergednodeid]=mergenodes(rn,connectivity,links,linksinconnect,fseg,mergenode1,mergenode2,MU,NU,a,Ec);    
                    %[links,connectivity,linksinconnect,fseg]=removelink(links,connectivity,linksinconnect,fseg,i);    
                end
            end
        elseif  (abs(rn(n1,lable))>surpos) & (abs( rn(n2,lable))>surpos )
            %mergenode1=n1;
            %mergenode2=n2;
            %[rn,connectivity,links,linksinconnect,fseg,mergednodeid]=mergenodes(rn,connectivity,links,linksinconnect,fseg,mergenode1,mergenode2,MU,NU,a,Ec);
            [links,connectivity,linksinconnect,fseg]=removelink(links,connectivity,linksinconnect,fseg,i);
            if (connectivity(n1,1)==0)
                deadnode=n1;
                [rn,connectivity,links]=removedeadnode(rn,connectivity,links,deadnode);
            elseif (connectivity(n2,1)==0)
                deadnode=n2;
                [rn,connectivity,links]=removedeadnode(rn,connectivity,links,deadnode);
            end
        elseif ( abs(rn(n1,lable))==surpos ) & ( abs(rn(n2,lable))==surpos ) 
            %if(connectivity(n1,1)<1)&(connectivity(n2,1)<1)
            %mergenode1=n1;
            %    mergenode2=n2;
            %    [rn,connectivity,links,linksinconnect,fseg,mergednodeid]=mergenodes(rn,connectivity,links,linksinconnect,fseg,mergenode1,mergenode2,MU,NU,a,Ec);
            %end
            [links,connectivity,linksinconnect,fseg]=removelink(links,connectivity,linksinconnect,fseg,i);
            if (connectivity(n1,1)==0)
                deadnode=n1;
                [rn,connectivity,links]=removedeadnode(rn,connectivity,links,deadnode);
            elseif (connectivity(n2,1)==0)
                deadnode=n2;
                [rn,connectivity,links]=removedeadnode(rn,connectivity,links,deadnode);
            end
        %   -----------------
        %       ._____.
        elseif ( abs(rn(n1,lable))<surpos ) & ( abs(rn(n2,lable))<surpos ) 
            if(connectivity(n1,1)==1) & (connectivity(n2,1)==1)
                if (rn(n1,lrn2)==0)| (rn(n2,lrn2)==0)
                    mergenode1=n1;
                    mergenode2=n2;
                    [rn,connectivity,links,linksinconnect,fseg,mergednodeid]=mergenodes(rn,connectivity,links,linksinconnect,fseg,mergenode1,mergenode2,MU,NU,a,Ec);    
                end        
            end    
        end
        
        i=i+1;
    end
            
           
function [links,connectivity,linksinconnect,fseg]=removelink(links,connectivity,linksinconnect,fseg,linkid);
% this subroutine is called by meshcoarsen
% this subroutine deletes the link information from the connectivity list
% removes the link from the link list and replaces the linkid with the last link
% after executing this subroutine all of the data structures should be clean

if(linkid>length(links(:,1)))
    disp(sprintf('try to remove link %d while total number of links is %d',linkid,length(links(:,1))));
    pause
end
% delete the linkid where it appears in connectivity list
% first appearance
nodeid1=links(linkid,1);
deadconnection1=linksinconnect(linkid,1);
[connectivity,linksinconnect]=removedeadconnection(connectivity,linksinconnect,nodeid1,deadconnection1);
% second appearance
nodeid2=links(linkid,2);
deadconnection2=linksinconnect(linkid,2);
[connectivity,linksinconnect]=removedeadconnection(connectivity,linksinconnect,nodeid2,deadconnection2);
% remove the link that no longer appears in the connectivity list and doesn't connect any nodes anymore
[links,connectivity,linksinconnect,fseg]=removedeadlink(links,connectivity,linksinconnect,fseg,linkid);



function [linksnew,connectivitynew,linksinconnectnew,fsegnew]=removedeadlink(linksnew,connectivitynew,linksinconnectnew,fsegnew,deadlink);
% this subroutine is called by meshcoarsen
% this subroutine replaces the link in linkid with the link in llinks
% repairs the connectivity and then deletes the llinks from links
% llinks should be the last linkid in links

% change the linkid in the connectivity to reflect the replacement
if(linksinconnectnew(deadlink,:)*linksinconnectnew(deadlink,:)'~=0)
    disp(sprintf('this dead link still has connections and should not be removed'));
    pause
end


llinks=length(linksnew(:,1));
if deadlink<llinks
    linksnew(deadlink,:)=linksnew(llinks,:);
    fsegnew(deadlink,:)=fsegnew(llinks,:);
    linksinconnectnew(deadlink,:)=linksinconnectnew(llinks,:);
    connectivitynew(linksnew(deadlink,1),2*linksinconnectnew(deadlink,1))=deadlink;
    connectivitynew(linksnew(deadlink,2),2*linksinconnectnew(deadlink,2))=deadlink;
end
linksnew(llinks,:)=[];
fsegnew(llinks,:)=[];
linksinconnectnew(llinks,:)=[];


function [rnnew,connectivitynew,linksnew]=removedeadnode(rnnew,connectivitynew,linksnew,deadnode);
% this subroutine is called by meshcoarsen and by removenode
% it removes nodes that are no longer part of the simulation
% and cleans up the data structures
% this subroutine replaces the node in i with the node in lrn
% repairs the links and then deletes the lrn from the node list
% lrn should be the last nodeid in rn
if(connectivitynew(deadnode,1)~=0)
    disp(sprintf('this deadnode still has connections and should not be removed'));
    deadnode
    rnnew
    connectivitynew
    pause
end

lrn=length(rnnew(:,1));
if deadnode<lrn
    rnnew(deadnode,:)=rnnew(lrn,:);
    connectivitynew(deadnode,:)=connectivitynew(lrn,:);
    for j=1:connectivitynew(deadnode,1) % change the nodeid in linksnew from lrn to i
        linksnew(connectivitynew(deadnode,2*j),connectivitynew(deadnode,2*j+1))=deadnode;
    end
end
rnnew(lrn,:)=[];
connectivitynew(lrn,:)=[];


function [connectivity,linksinconnect]=removedeadconnection(connectivity,linksinconnect,nodeid,deadconnection);
%This subroutine deletes an entry in a node's connectivity list and updates the linksinconnet array

lastconnection=connectivity(nodeid,1);
%remove the entry in linksinconnect to show that the connectivity data no longer exits for that link
linksinconnect(connectivity(nodeid,2*deadconnection),connectivity(nodeid,2*deadconnection+1))=0;
if (lastconnection>deadconnection)
   %replace link in the connectivitylist with the lastlink in the connectivity list
   connectivity(nodeid,2*deadconnection:2*deadconnection+1)=connectivity(nodeid,2*lastconnection:2*lastconnection+1);
   % update linksinconnect to reflect the change in position of the lastconnection
   linksinconnect(connectivity(nodeid,2*deadconnection),connectivity(nodeid,2*deadconnection+1))=deadconnection;
end
connectivity(nodeid,2*lastconnection:2*lastconnection+1)=[0 0];
connectivity(nodeid,1)=lastconnection-1;