function [miu,miu_c,Kvc,KEC,KC,gammaMiu_el] = FEMcoupler_new(rn,links,fn_c,normc,xnodes,mno,wx,wy,wz,mx,my,mz,mel,nc,ke,Dv,nodelist,conlist,connectivity)
    
%26/11/2020 By fengxian Liu
% solving the bulk diffusion equation with consideration of the core region

%----------- get gammaU (the set of nodes on dislocation core) and U_dis (potential on dislocation core) from the location of dislocation lines rn, links
[gammaMiu_num, gammaMiu_el, gammaMiu_area] = findcore_new(rn,links,xnodes,wx,wy,wz,mx,my,mz,mel,nc);
%gammaMiu_num: the number of element with each segment
%gammaMiu_el: element number of each segment
%gammaMiu_area: integration of s inside of each element

%--------------- calculate ke in the core region
nseg=size(links,1);
nnode=size(rn,1);

Ke=zeros(8,8,mel);

Kc=zeros(2,2,nseg);
kvc=zeros(2,2,nseg);

KEC=sparse(mno,nnode);
% KEC1=sparse(mno,nnode);

r0_rc=50; %r0/rc
for i=1:nseg
    n0=links(i,1);
    n1=links(i,2);
    rt=rn(n1,1:3)-rn(n0,1:3);
    l=norm(rt);
    
    %calculate kvc
    burgs=links(i,3:5);    
    if l>0        
        theta=dot(rt,burgs)/norm(burgs)/l;
        theta=acos(theta);
        sintheta=abs(sin(theta));
    else
        sintheta=0;
    end
    kvc0=[1/3 1/6;1/6 1/3];
    if sintheta < sin(5*pi/180) %  cutt off angle 5;
        kvc(:,:,i)=kvc0*1e5*l;
    else
        kvc(:,:,i)=kvc0*l*sintheta;
    end
 
    %   
    for n=1:gammaMiu_num(i)
        s0=gammaMiu_area(i,2*n-1);
        s1=gammaMiu_area(i,2*n);
        element=gammaMiu_el(i,n);
        
        %calculate the elemental matrix for the near core element 
        [kee,kec,kcc]=intNs_gauss(rn(n0,1:3),rn(n1,1:3),xnodes(nc(element,1),1:3),s0,s1,wx,wy,wz);
        
        alpha0=pi*Dv*l/log(r0_rc);
        kee=alpha0*kee;
        kec=alpha0*kec;
        kcc=alpha0*kcc;
        
        Ke(:,:,element) =Ke(:,:,element) + kee;
       
        Kc(:,:,i)=Kc(:,:,i)+ kcc;
            
        dof=nc(element,1:8);
        
        KEC(dof,n0)=KEC(dof,n0)+kec(:,1);
        KEC(dof,n1)=KEC(dof,n1)+kec(:,2);
        
    end
end

% disp('global stifness matrix assembly ...')
%KE, mno*mno, golbal matrix for all the FE nodes
%KC, nnode*nnode, gobal matrix for all the dislocation nodes
%KEC, mno*nnode, complmentary matrix to relate miu and miu_c

% KE=zeros(mno,mno);
% KC=zeros(nnode,mno);
%  KEC=zeros(mno,nnode);

a=1:8; %local node numbers 
dofLocal=a;
% ntriplets = mno;
I = zeros(mno,1);
J = zeros(mno,1);
X = zeros(mno,1);
ntriplets = 0;

%matrix KE
elements=unique(gammaMiu_el);
for p =1:mel
    gn=nc(p,a); % global node numbers
    dof(1:8)=gn; % global degree of freedom
    for i =1:8
        for j =1:8
            ntriplets = ntriplets + 1;
            I(ntriplets) = dof(i);
            J(ntriplets) = dof(j);
            if any(elements==p)
                X(ntriplets) = Ke(dofLocal(i),dofLocal(j),p);
            else
                X(ntriplets) = ke(dofLocal(i),dofLocal(j),p);
            end
        end
    end
end
KE= sparse(I,J,X,mno,mno); %the (full) global stiffness matrix for FE nodes


%matrix KC
ntriplets=0;
I = zeros(nnode,1);
J = zeros(nnode,1);
X = zeros(nnode,1);
Y = zeros(nnode,1);
a=1:8; %local node numbers 
dofLocal=a;
for k=1:nseg
    n0=links(k,1);
    n1=links(k,2);
    dof1(1:2)=[n0 n1];
    for i=1:2
        for j=1:2
          ntriplets=ntriplets+1;
          I(ntriplets) = dof1(i);
          J(ntriplets) = dof1(j);
          X(ntriplets) = Kc(dofLocal(i),dofLocal(j),k);
          Y(ntriplets) = kvc(dofLocal(i),dofLocal(j),k);
        end        
            
    end
end
KC= sparse(I,J,X,nnode,nnode); %the (full) global stiffness matrix for dislocation nodes
Kvc= sparse(I,J,Y,nnode,nnode); 


%{miu_c}
miu_c=zeros(nnode,1);
for i=1:nnode
    miu_c(i)=norm(fn_c(i,1:3))/1; %miu_c=fc/b, N/m^2
    if dot(normc(i,:),fn_c(i,:))<0
         miu_c(i)= -miu_c(i);
    end
end


%--------------------------- solve function 
%[KE]{miu}+[KEC]{miu_c}={0}  [KE]: mno*mno, m^3*mol/J*s ; {miu}: mno*1, J/m^3;  [KEC]:mno*nnode; {miu_c}: nnode*1; {0}: mno*1, matter, mol/s 
miu=zeros(mno,1);
J=zeros(mno,1); 
J = J-KEC*miu_c; %
miu = KE\J; %J/m^3


%--------------------------- calculate {vc}
%[KEC]'*{miu}+[KC]*{miu_c}=[Kvc]*{vc}
%[Kvc]
end