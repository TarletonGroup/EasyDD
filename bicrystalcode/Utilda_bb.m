function [Ux, Uy, Uz] = Utilda_bb(rn,links,gnl,nu,xnodes,dx,dy,dz,mx,my,mz)

%This function finds the displacements on the FEM nodes in node list gnl
%caused by all dislocation segments in links. Can be adapted to use
%segments and x0 in place of rn, links, gnl and xnodes. However the
%correction for erronious nodes has not been made, and segments cannot be
%used in place of rn and links once it is

nodenum=size(gnl,1);
segnum=size(links,1);
Utilda=zeros(nodenum,3);
x=dx/mx;
y=dy/my;
z=dz/mz;
C=[x;y;z];

for i=1:nodenum
   nodepoint=xnodes((gnl(i)),1:3);
   p=nodepoint';
   for j=1:segnum
      A=rn(links(j,1),1:3)';
      B=rn(links(j,2),1:3)';
      b=links(j,3:5)';
      Utilda(i,:)=Utilda(i,:)+displacement_et(p,A,B,C,b,nu);
       
   end
   
end

Ux=Utilda(:,1);
Uy=Utilda(:,2);
Uz=Utilda(:,3);

end