function [ep_inc,wp_inc]=calcplasticstrainincrement_HY(rnnew,rn,links,Volume)
seg=      rn(links(:,2),1:3)   - rn(links(:,1),1:3);
segnew=rnnew(links(:,2),1:3)- rnnew(links(:,1),1:3);

%HY20192011: remove the virtual nodes
virtual = [find(rn(links(:,1),end)==67);find(rn(links(:,2),end)==67);...
    find(rnnew(links(:,1),end)==67);find(rnnew(links(:,2),end)==67)];
virtual = unique(virtual,'stable');  
% seg(virtual,:) = 0;
% segnew(virtual,:) = 0;

dx1=rnnew(links(:,2),1:3)-rn(links(:,1),1:3);
dA=cross(segnew+seg,dx1);
dA(virtual,:) = 0;
fp_inc=0.5.*(links(:,3:5)'*dA)./Volume;
ep_inc=0.5.*(fp_inc+fp_inc');
wp_inc=0.5.*(fp_inc-fp_inc');

%HY20191103: corrections
% fp_inc = zeros(3);
% for i=1:size(links,1)
%     fp_inc = fp_inc+(links(i,3:5)'*links(i,6:8)*dA(i));
% end
% fp_inc=0.5.*fp_inc./Volume;
% ep_inc=(fp_inc+fp_inc');
% wp_inc=0.5.*(fp_inc-fp_inc');