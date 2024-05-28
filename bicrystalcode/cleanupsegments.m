function [rnnew,linksnew,connectivitynew,linksinconnectnew,fsegnew] = cleanupsegments(rn,links,fseg)

%This functions removes segments with effective zero Burgers vectors and
%any floating nodes generated as a result.

bvecs=links(:,3:5);
bvecs(abs(bvecs)<1E-2)=0;
links(:,3:5)=bvecs;
fseg(~any(links(:,3:5),2),:)=[];
links(~any(links(:,3:5),2),:)=[];

[connectivity,linksinconnect,~]=genconnectivity(rn,links,4);

if any(connectivity(:,1)==0)
    for i=1:size(rn,1)
        if connectivity(i,1)==0
            links(links(:,1)>i,1)=links(links(:,1)>i,1)-1;
            links(links(:,2)>i,2)=links(links(:,2)>i,2)-1;
        end
    end
end
rn(connectivity(:,1)==0,:)=[];
connectivity(connectivity(:,1)==0,:)=[];


rnnew=rn;
linksnew=links;
connectivitynew=connectivity;
linksinconnectnew=linksinconnect;
fsegnew=fseg;
