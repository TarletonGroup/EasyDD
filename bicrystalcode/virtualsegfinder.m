%This script runs through the rn and links arrays and finds the segments
%that consist of a surface node (flagged 6) and a virtual node (flagged 67)
%We need these segments to correct for the stress at the surface, since we
%ignore all other virtual segments (thereby not having a closed integral).
%This is acceptable as long as the surf-virt nodes are 100x-1000x longer
%than the dislocation segment exiting the surface.

%virtID index says whether 1st or 2nd entry in links is surface node.
%realID index says whether 1st or 2nd entry in links is surface node.

%realsurflist gives links real links to surface nodes. This is useful to
%know when correcting for forces at surface nodes.

%function [surfvirtlist,virtIDlist,realsurflist,realIDlist] = virtualsegfinder(rn,links)

function [index,index2,indexR,indexR2,segpairs] = virtualsegfinder(rn,links,segments)

Slinks = size(links,1);

index_a=(rn(links(:,1),end)==6 & rn(links(:,2),end)==67);  %finds links attaching surface nodes to vitual nodes
index_b=(rn(links(:,2),end)==6 & rn(links(:,1),end)==67);
index = (index_a | index_b);                               %index for said links
index2 = zeros(Slinks,1);                                  %flags whether surface node is 1st or 2nd in links
index2(index_a==true) = 1;
index2(index_b==true) = 2;
indexR_a=(rn(links(:,1),end)==6 & rn(links(:,2),end)==0);  %finds links attaching surface nodes to real nodes
indexR_b=(rn(links(:,2),end)==6 & rn(links(:,1),end)==0);
indexR_c=(rn(links(:,1),end)==6 & rn(links(:,2),end)==7);
indexR_d=(rn(links(:,2),end)==6 & rn(links(:,1),end)==7);
indexR = (indexR_a | indexR_b | indexR_c | indexR_d);                            %index for said links
indexR2 = zeros(Slinks,1);                                 %flags whether surface node is 1st or 2nd in links
indexR2(indexR_a==true | indexR_c==true) = 1;
indexR2(indexR_b==true | indexR_d==true) = 2;

index2(index2==0) = [];                                    %removes zeros from index2 and indexR2
indexR2(indexR2==0) = [];

virtsegs=segments(index,:);                                %lists segments attached to surface nodes that are virtual
realsegs=segments(indexR,:);                               %lists segments attached to surface nodes that are real
SR = size(realsegs,1);                                     %number of real segments attached to surface nodes
SV = size(virtsegs,1);                                     %number of virtual segments attached to surface nodes

segpairs=zeros(2*SV*SR,14);                                %preallocates memory, accounting for multiple real and/or virtual segments attached to one surface node

for i=1:SR
    for j=1:SV
        if virtsegs(j,index2(j)) == realsegs(i,indexR2(i)) %virtual seg is paired with real seg
            segpairs(2*(i-1)*SV+2*j-1,:) = realsegs(i,:);  %first seg in pair is real
            segpairs(2*(i-1)*SV+2*j,:) = virtsegs(j,:);    %second seg in pair is virtual
        end
    end
end

segpairs(~any(segpairs,2),:) = [];                         %removes empty rows from segpairs

end
