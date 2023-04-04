hold on
for n=1:L1
%      quiver3(rn(n,1),rn(n,2),rn(n,3),normc(n,1),normc(n,2),normc(n,3),1e2,'-b','LineWidth',1);
      quiver3(rn(n,1),rn(n,2),rn(n,3),vn_b(n,1),vn_b(n,2),vn_b(n,3),2*1e12,'r','LineWidth',1);
%       quiver3(rn(n,1),rn(n,2),rn(n,3),vn_bulk(n,1),vn_bulk(n,2),vn_bulk(n,3),1e2,'r','LineWidth',1);
%      quiver3(rn(n,1),rn(n,2),rn(n,3),vn_core(n,1),vn_core(n,2),vn_core(n,3),1e2,'b','LineWidth',1);
%      quiver3(rn(n,1),rn(n,2),rn(n,3),fn_c(n,1),fn_c(n,2),fn_c(n,3),3e-5,'c','LineWidth',1);
    hold on
end


hold on
for n=1:nnode
    quiver3(rn(n,1),rn(n,2),rn(n,3),miu_c(n,1),3e-7,'g','LineWidth',1);

    hold on
end