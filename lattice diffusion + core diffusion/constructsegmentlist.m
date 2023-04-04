function segments=constructsegmentlist(rn,links)
% 
[LINKMAX,~]=size(links);

segments=zeros(LINKMAX,14);

if sum(any(links(:,1:2)==0)==1)~=0
%    tic 
    nseg=0;
    for i=1:LINKMAX
        n0=links(i,1);
        n1=links(i,2);
        if((n0~=0)&&(n1~=0))
            nseg=nseg+1;
            segments(nseg,:)=[links(i,1:5),rn(n0,1:3),rn(n1,1:3),links(i,6:8)];
        end
    end
    
    segments=segments(1:nseg,:);
%     toc
else
%     tic
   segments(:,1:5)=links(:,1:5);
   segments(:,6:8)=rn(links(:,1),1:3);
   segments(:,9:11)=rn(links(:,2),1:3);
   segments(:,12:14)=links(:,6:8);
%    toc
end
% [LINKMAX,~]=size(links);
% 
% segments=zeros(LINKMAX,14);
% nseg=0;
% for i=1:LINKMAX
%     n0=links(i,1);
%     n1=links(i,2);
%     if(     (        n0     ~= 0  &&    n1     ~= 0  ))
%           if (rn(n0, end) == 67 || rn(n1, end) == 67 )
%               continue
%           end
%         nseg=nseg+1;
%         segments(nseg,:)=[links(i,1:5),rn(n0,1:3),rn(n1,1:3),links(i,6:8)];
%     end
% end
% segments=segments(1:nseg,:);