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

function [index,index2,indexR,indexR2,segpairs] = virtualsegfinder_old(rn,links)

surfnodes = find(rn(:,4)==6); %finds surface nodes

S = length(surfnodes);
Slinks = size(links,1);

index = false(Slinks,1);
index2 = zeros(Slinks,1);
indexR = false(Slinks,1);
indexR2 = zeros(Slinks,1);

for i = 1:S
    [row,col] = find(links(:,1:2) == surfnodes(i));
    S_conns = length(row);

    %for each connection to surface node, check whether virtual or real
    for j = 1:S_conns
        if col(j) == 1
            link = links(row(j),2); %2nd row
            if rn(link,4) == 67 %virtual node flag
                index(row(j)) = true;
                index2(row(j)) = 1; %surf node 1st col of links array
            elseif rn(link,4) == 0 %real node flag
                indexR(row(j)) = true;
                indexR2(row(j)) = 1; %surf node 1st col of links array
            end
        else
            link = links(row(j),1); %1st row
            if rn(link,4) == 67 %virtual node flag
                index(row(j)) = true;
                index2(row(j)) = 2; %surf node 2nd col of links array
            elseif rn(link,4) == 0 %real node flag
                indexR(row(j)) = true;
                indexR2(row(j)) = 2; %surf node 1st col of links array
            end
        end
    end

end

index2(index2==0) = [];
indexR2(indexR2==0) = [];

segments=constructsegmentlist(rn,links);
segpairs=segpairing(segments,index,indexR,indexR2);

segpairs(segpairs==0)=[];

end

function segpairs=segpairing(segments,index,indexR,indexR2)

virtsegs=segments(index,:);
realsegs=segments(indexR,:);

SR = size(realsegs,1);
SV = size(virtsegs,1);

segpairs=zeros(SV*SR,14); %first is real, second is virtual

counter=1;
for i=1:SR
    for j=1:SV
    condition = sum(virtsegs(j,1:2) == realsegs(i,indexR2(i)));
        if condition == 1 %virtual seg is paired with real seg
            segpairs(2*counter-1,:) = realsegs(i,:);
            segpairs(2*counter,:) = virtsegs(j,:);
            counter=counter+1;
        end
    end
end
end
