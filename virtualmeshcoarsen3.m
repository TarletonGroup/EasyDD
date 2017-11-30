%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check whether virtual nodes can be eliminated based on:
% 1) If they are not connected to any surface nodes
% 2) If they are not due to an angle change in the simulated volume surface
% 
% Bruce Bromage
% Michromechanical Testing Group
% Department of Materials, University of Oxford
% bruce.bromage@materials.ox.ac.uk
% May 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew]=virtualmeshcoarsen3(rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew,MU,NU,a,Ec)

%This function can coarsen away cross slip at the surface and may need to be
%improved

lcrit=1e3;   %this is an arbitrary distance and should be corrected to be related to rmax
acrit=0.5*(lcrit^2)*sin(2*pi/360);   %this is an area related to a desired resolution of angle change and can be altered as required using the sine term
node_id=1;

while node_id<=size(rnnew,1)    %start from the top of rn and go through all nodes
    
    if connectivitynew(node_id,1)==2 && rnnew(node_id,end)==67   %target only virtual nodes with two connections
            link_node1=linksnew(connectivitynew(node_id,2),3-connectivitynew(node_id,3));   %first node linked to target node
            link_node2=linksnew(connectivitynew(node_id,4),3-connectivitynew(node_id,5));   %second node linked to target node
            
            if rnnew(link_node1,end)==67 && rnnew(link_node2,end)==67   %only if both nodes linked to target node are virtual
                link_vec1=rnnew(link_node1,1:3)-rnnew(node_id,1:3);   %vector of link 1
                link_vec2=rnnew(link_node2,1:3)-rnnew(node_id,1:3);   %vector of link 2
                angle=acos(norm(dot(link_vec1,link_vec2))/(norm(link_vec1)*norm(link_vec2)));   %angle between link 1 and link 2
                area=0.5*norm(link_vec1)*norm(link_vec2)*sin(angle);   %area of triangle formed by link 1 and link 2
                
                if norm(link_vec1)<lcrit && area<acrit   %if length of link 1 and the angle change are below critical size then merge
                    [rnnew,connectivitynew,linksnew,linksinconnectnew,fsegnew,~]=mergenodes(rnnew,connectivitynew,linksnew,linksinconnectnew,fsegnew,link_node1,node_id,MU,NU,a,Ec);
                    
                elseif norm(link_vec2)<lcrit && area<acrit  %if length of link 2 and the angle change are below critical size then merge
                    
                    [rnnew,connectivitynew,linksnew,linksinconnectnew,fsegnew,~]=mergenodes(rnnew,connectivitynew,linksnew,linksinconnectnew,fsegnew,link_node2,node_id,MU,NU,a,Ec);
                    
                else
                   node_id=node_id+1;   %target next node if both links are larger than critical size
                end
                
            else
                node_id=node_id+1;   %target next node if at least one of the links is not virtual
            end
            
    else
        node_id=node_id+1;   %target next node if current target node is not virtual or does not have two connections
    end
    
end

end
            