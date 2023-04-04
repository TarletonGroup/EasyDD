function [vnvec,fseg] = drndt(rnvec,flag,MU,NU,a,Ec,links,connectivity,appliedstress,...
                 mobility,SimBox,xnodes,mno,wx,wy,wz,mx,my,mz,mel,nc,ke,Dv)
                              

%unscramble rn
rn = [reshape(rnvec,length(rnvec)/3,3),flag];



%nodal driving force
fseg = segforcevec(MU,NU,a,Ec,rn,links,appliedstress,0);  
%mobility function
% fseg
%fseg = imageforce(MU,NU,a,rn,links,connectivity,fseg,Surpos,Btype);

%mobility function
[vn,vn_c,fn] = feval(mobility,fseg,rn,links,connectivity,xnodes,mno,wx,wy,wz,mx,my,mz,mel,nc,ke,Dv,[],[]);

% [vn1,fn_g] = feval(mobility_g,fseg,rn,links,connectivity,[],[])
%fixed nodes (flag~=0) are not allowed to move
% vn=vn+vn1;
vn=vn.*((rn(:,4)==0)*[1 1 1]);






%make up the last
vnvec=reshape(vn,length(vn(:,1))*3,1);

