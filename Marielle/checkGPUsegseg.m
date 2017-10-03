

figure(1);
clf
%amag = 3.18e-4; %lattice vector BCC W, in microns
amag=1;
plot3(0,0,0); hold on;
LINKMAX=length(links(:,1));
for i=1:LINKMAX,
    n0=links(i,1);
    n1=links(i,2);
    %to skip external nodes...
    if rn(n0,end)==67||rn(n1,end)==67
       continue;
    end
    lvec = amag*(rn(n1,1:3)-rn(n0,1:3));
    plane_n = links(i,6:8);
    bvec = links(i,3:5);
    plane_n = cross(lvec/norm(lvec),bvec/norm(bvec));
%     if plane_n == [-1 0 1]/sqrt(2)
    r0 = rn(n0,1:3)*amag;
%      if rn(n0,4)==0
        %filter out "infinity" lines
       plot3(rn([n0,n1],1)*amag,rn([n0,n1],2)*amag,rn([n0,n1],3)*amag,'r','LineWidth',2);  
        quiver3(r0(1),r0(2),r0(3),lvec(1),lvec(2),lvec(3),'r','LineWidth',1); 
        %plot3(rn(n0,1)*bmag,rn(n0,2)*bmag,rn(n0,3)*bmag,'k.');            
        %plot3(rn(n1,1)*bmag,rn(n1,2)*bmag,rn(n1,3)*bmag,'k.');            
%      end 
  
end
 hold on; 
%amag = 3.18e-4; %lattice vector BCC W, in microns
amag=1;
plot3(0,0,0); hold on;
LINKMAX=length(l(:,1));
for i=1:LINKMAX,
    n0=l(i,1);
    n1=l(i,2);
    %to skip external nodes...
    if r(n0,end)==67||r(n1,end)==67
       continue;
    end
    lvec = amag*(r(n1,1:3)-r(n0,1:3));
    plane_n = l(i,6:8);
    bvec = l(i,3:5);
    plane_n = cross(lvec/norm(lvec),bvec/norm(bvec));
%     if plane_n == [-1 0 1]/sqrt(2)
    r0 = r(n0,1:3)*amag;
%      if rn(n0,4)==0
        %filter out "infinity" lines
       plot3(r([n0,n1],1)*amag,r([n0,n1],2)*amag,r([n0,n1],3)*amag,'b','LineWidth',2);  
        quiver3(r0(1),r0(2),r0(3),lvec(1),lvec(2),lvec(3),'r','LineWidth',1); 
        %plot3(rn(n0,1)*bmag,rn(n0,2)*bmag,rn(n0,3)*bmag,'k.');            
        %plot3(rn(n1,1)*bmag,rn(n1,2)*bmag,rn(n1,3)*bmag,'k.');            
%      end 
  
end
