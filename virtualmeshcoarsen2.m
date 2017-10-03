%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check whether virtual nodes can be eliminated based on:
% 1) Segment length (i.e. between lmin lmax)
% Simpler/faster than virtualmeshcoarsen.m
%
% Francesco Ferroni
% Defects Group / Materials for Fusion and Fission Power Group
% Department of Materials, University of Oxford
% francesco.ferroni@materials.ox.ac.uk
% February 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [rnnew,linksnew,connectivitynew,linksinconnectnew] = virtualmeshcoarsen2(rn,links,maxconnections,lmin)

%| rn(loop_list(:,1),4)==6
% N.B. loop_list gives a list of nodes that below in a loop, in sequence.
% However, it does not give the "direction" of this sequence since segments
% are vectors, hence you need to check this based on input links array.

lrn=size(rn,1); %total number of nodes M 
[connectivity,linksinconnect,~]=genconnectivity(rn,links,maxconnections);
counter=1;
fseg=[]; %empty

while counter>0;
    counter=0;
    i=1;
    while i<lrn % i<total numbers of nodes M
        if rn(i,end)~=67 || connectivity(i,1)~=2
            %ignore domain/fixed/surface nodes
            %ignore junctions
            i=i+1;
            continue;
        end
        node1 = links(connectivity(i,2) , logicswap(connectivity(i,3)));
        node2 = links(connectivity(i,4) , logicswap(connectivity(i,5)));
        seg1 = rn(i,1:3) - rn(node1,1:3);
        seg2 = rn(i,1:3) - rn(node2,1:3);
        length = sqrt(sum(seg1.*seg1)) + sqrt(sum(seg2.*seg2)); %length=seg1+seg2 ? M

        if length < lmin
            %check angle
            %angle=acos(dot(seg1,seg2)/norm(seg1)/norm(seg2));
            %if angle < (pi/10) || angle > (9*pi/10) %if the curvature is very high, then don't delete.
            % remove node and remesh
            counter=counter+1;
            [rn,connectivity,links,linksinconnect,~,~]=...
                mergenodes(rn,connectivity,links,linksinconnect,fseg,node1,i);
            %end
        end
        lrn=size(rn,1);
        i=i+1;
    end
    
    %clean-up orphan nodes flagged as superfluous and re-index nodal structure.
    %[rn,links] = clearorphannodes(rn,links);
    %update connectivity and linksinconnect.
    %[~,~,~] = genconnectivity(rn,links,maxconnections);  
end

%clean-up orphan nodes flagged as superfluous and re-index nodal structure.
[rnnew,linksnew] = clearorphannodes(rn,links);
%update connectivity and linksinconnect.
[connectivitynew,linksinconnectnew] = genconnectivity(rn,links,maxconnections);

end

function [rnnew,linksnew] = clearorphannodes(rnnew,linksnew)

orphan_nodes = rnnew(:,end)==-1;
while any(orphan_nodes)
    for i=1:size(rnnew,1)    
        if rnnew(i,4)==-1
            [row,col] = find(linksnew(:,1:2)==i);  
            %check it's a pair!
            if length(row)~=2
                disp('WARNING! Not analyzing pairs! See virtualmeshcoarsen.m (row 79)')
            end 
            if col(1 )==2
                %merge
                linksnew(row(1),2) = linksnew(row(2),2); 
                linksnew(row(2),:) = [];
            elseif col(1)==1
                %merge
                linksnew(row(1),1) = linksnew(row(2),1); 
                linksnew(row(2),:) = [];
            end
            rnnew(i,end)=1; %change flag
            break;
        end    
    end
%update orphan_nodes list
orphan_nodes = (rnnew(:,end)==-1);
end

orphan_nodes = (rnnew(:,end)==1); %look for flagged rn.
index_corr = cumsum(orphan_nodes);
linksnew(:,1) = linksnew(:,1) - index_corr(linksnew(:,1)); %update numbering
linksnew(:,2) = linksnew(:,2) - index_corr(linksnew(:,2)); %update numbering
rnnew(orphan_nodes,:) = []; %finally get rid of nodes

end

function out = logicswap(in)
if in==1
    out=2;
elseif in==2;
    out=1;
else
    disp('Incorrect use of logicswap function');
end
end