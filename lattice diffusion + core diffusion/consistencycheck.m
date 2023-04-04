function c=consistencycheck(rn,links,connectivity,linksinconnect,Box)
%self consistency check for topology

if(max(max(connectivity))~=length(links(:,1)))
    disp(sprintf('max connectivity = %d length(links) = %d',max(max(connectivity)),length(links(:,1))))
    pause
end

nnodes=length(rn(:,4));
for i=1:nnodes,
    numNbrs=connectivity(i,1);
    if(sum(connectivity(i,2*(numNbrs+1):end)~=0))
        disp('inconsistent connectivity list');
        i
        connectivity(i,:)
        pause
    end
end
% check for conservation of burgers vector at nodes
%modified: check only the node are in the box
for i=1:nnodes
    if (rn(i,4)~=7)&(abs(rn(i,1))<Box(1))&(abs(rn(i,2))<Box(2))&(abs(rn(i,3))<Box(3))
        totalb=zeros(1,3);
        numNbrs=connectivity(i,1);
        for j=1:numNbrs
        linkid=connectivity(i,2*j);
        posi=connectivity(i,2*j+1);
        totalb=totalb+(3-2*posi).*links(linkid,3:5);  
        end
        if totalb*totalb'~=0
            disp(sprintf('the Burgers vector is not conserved at node %d',i));
            %totalb*totalb'
            pause       %some segment may be at the surface
        end
    end
end

%if(sum(sum(connectivity==length(links(:,1))))~=2)
%    disp(sprintf('arm %d has %d ends',length(links(:,1)),sum(sum(connectivity==length(links(:,1))))))
%    pause
%end

for i=1:nnodes,
    numNbrs=connectivity(i,1);
    for j=1:numNbrs,
        for k=j+1:numNbrs,
            if(connectivity(i,2*j)==connectivity(i,2*k))
                connectivity(i,:)
                disp(sprintf('node %d connected with link %d twice',i,connectivity(i,2*j)));
                pause
            end
        end
    end
end

for i=1:nnodes,
    numNbrs=connectivity(i,1);
    for j=1:numNbrs,
        n0=links(connectivity(i,2*j),connectivity(i,2*j+1));    
        if(i~=n0)
            disp(sprintf('connectivity(%d,:)',i));
            connectivity(i,:)
            disp(sprintf('links(%d,:)=(%d,%d)',connectivity(i,2*j),links(connectivity(i,2*j),1),links(connectivity(i,2*j),2)));
        end
    end
end

for i=1:length(links(:,1))    
    j=connectivity(links(i,1),2*linksinconnect(i,1));
    k=connectivity(links(i,2),2*linksinconnect(i,2));
    if(i~=j)|(i~=k)
        disp(sprintf('inconsistent link %d %d %d',i,j,k));
        disp(sprintf('links(%d,:)=(%d,%d)',i,links(i,1),links(i,2)));
        disp(sprintf('linksinconnect(%d,:)=(%d,%d)',i,linksinconnect(i,1),linksinconnect(i,2)));
        disp(sprintf('connectivity(%d,:)=(%d %d %d %d %d %d %d %d %d)',links(i,1),...
            connectivity(links(i,1),1),connectivity(links(i,1),2),connectivity(links(i,1),3),...
            connectivity(links(i,1),4),connectivity(links(i,1),5),connectivity(links(i,1),6),...
            connectivity(links(i,1),7),connectivity(links(i,1),8),connectivity(links(i,1),9) ));
        disp(sprintf('connectivity(%d,:)=(%d %d %d %d %d %d %d %d %d)',links(i,2),...
            connectivity(links(i,2),1),connectivity(links(i,2),2),connectivity(links(i,2),3),...
            connectivity(links(i,2),4),connectivity(links(i,2),5),connectivity(links(i,2),6),...
            connectivity(links(i,2),7),connectivity(links(i,2),8),connectivity(links(i,2),9) ));
            
        pause
    end
end
