function [connectivity,linksinconnect,node_connectivity]=genconnectivity(rn,links,maxconnections)
%the connectivity is basically the inverse of the links array
% for each node it the connectivity tells how many links that node is part of 
% and then lists the link indices of the system
% connectivity(i,1)= number of links that node i is part of 
% connectivity(i,2:7)= the list of connectivity(i,1) links that node i is part of
rnlength=length(rn(:,4));
connectivity=zeros(rnlength,1+2*maxconnections);
linkslength=length(links(:,1));
linksinconnect=zeros(linkslength,2);  %linksinconnect has 2 columns, one for node 1 of linkid one for node 2 of linkid M
for i=1:linkslength
    if links(i,1)~=0 %node 1 different from 0
        a=links(i,1); %node1 M 
        b=links(i,2); %node2 M 
        connectivity(a,1)=connectivity(a,1)+1;
        connectivity(b,1)=connectivity(b,1)+1;
        connectivity(a,2*connectivity(a,1):2*connectivity(a,1)+1)=[i 1]; %[i 1] i=linkid , 1=first node of the linkid M
        connectivity(b,2*connectivity(b,1):2*connectivity(b,1)+1)=[i 2]; %[i 2] i=linkid , 2=second node of the linkid M
        linksinconnect(i,1)=connectivity(a,1); %linksinconnect=connectivity of the node
        linksinconnect(i,2)=connectivity(b,1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_nodes = size(rn,1);
node_connectivity = zeros(num_nodes,maxconnections);
%node_connectivity : each line = 1 node, the numbers correspond to the
%other nodes the node in question is linked
for i=1:num_nodes
        for j=1:(size(connectivity,2)-1)*0.5   %(size(connectivity,2)-1)*0.5=maxconnection...M
            if connectivity(i,2*j) == 0
                continue;
            else
            pos=connectivity(i,2*j+1);
            if pos==1
                pos=2;
            else 
                pos=1;
            end
            
            node_connectivity(i,j) = links(connectivity(i,2*j),pos);
            
            end
        end
end

end