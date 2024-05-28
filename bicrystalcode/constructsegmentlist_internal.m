function segments=constructsegmentlist_internal(rn,links)

[LINKMAX,~]=size(links);

segments=zeros(LINKMAX,14);
nseg=0;
for i=1:LINKMAX
    n0=links(i,1);
    n1=links(i,2);
    if((n0~=0)&&(n1~=0)&&(rn(n0,end)~=67)&&(rn(n1,end)~=67))
        nseg=nseg+1;
        segments(nseg,:)=[links(i,1:5),rn(n0,1:3),rn(n1,1:3),links(i,6:8)];
    end
end
segments=segments(1:nseg,:);