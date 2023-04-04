% modified by gaoyuan 09-09-18

function [v_c,fn_c] = mobbcc_core(fseg,rn,links,connectivity,nodelist1,conlist1)

global bur tem

R = 8.314; % J/mol/K
Boltz    = 1.3807e-23; % J/K
Atomvol  = (2.856)^3*1e-30;  %volum of the atom, m^3/s
rc       = 5*bur*1e-9;  % dislocation core radius, m

%%-----------------------------------------------
pre_exponential_c = 1e-23; % ac*Dc0  m^4/s
active_energy = 174*1e3; % kJ/mol  1.8034eV
D_core = pre_exponential_c*exp(-active_energy/R/tem);  % core diffusivity
Dc=D_core*Atomvol/Boltz/tem;
 
%----------------------------------------------------------------------------------------
L1 = size(rn,1);    % when mobility law is call form remesh fucnction, L1~=L0;
nodelist = linspace(1,L1,L1)';
[L2,L3] = size(connectivity);
conlist = zeros(L2,(L3-1)/2+1);
conlist(:,1) = connectivity(:,1);
for i=1:L2
    connumb = conlist(i,1);
    conlist(i,2:connumb+1) = linspace(1,connumb,connumb);
end

%% dislocation climb dominated by core diffusion  ( edited by Ally  12/11/2018 )

v_c=zeros(L1,3);  % nodal climb velocity
fn_c=zeros(L1,3); % nodal climb force
lamla_c=zeros(L1,1); %chemical potential at each node

%%------------------------------------------------------------ node force
nseg=size(links,1);
    nnode=size(rn,1);
    C=zeros(nnode,nnode+nseg);

    k=zeros(nseg+nnode,nseg+nnode);
    fnc=zeros(nnode,1);
    fn_c=zeros(nnode,3);
    K=zeros(nnode*2+nseg,nnode*2+nseg);
    l=zeros(nseg,3); 
    nc=zeros(nnode,3);
    fn=zeros(nnode,3);
    Fnc=zeros(nnode+nseg,1);
    F=zeros(nnode*2+nseg);
    
    linvec=zeros(nnode,3);
    linel=zeros(nnode,1);
    vt=zeros(nnode,3); % line velocity;

    
    cc=[-3/8,-1/8, -1;-1/8,-3/8, 1];
    %%%%%%calculate the climb direction and climb force for each node 
    L=0;
    for n=1:L1
        n0=nodelist(n);
        numNbrs=conlist(n,1);
        rt=zeros(numNbrs,3);
        
        for i=1:numNbrs
            ii=conlist(n,i+1);
            linkid=connectivity(n0,2*ii);
            posinlink=connectivity(n0,2*ii+1);
            n1=links(linkid,3-posinlink);
            rt=rn(n1,1:3)-rn(n0,1:3);
            
            
            linvec(n,:)= linvec(n,:)+rt*(-1)^(3-posinlink);
            l(linkid)=norm(rt);
            linel(n)=linel(n)+l(linkid)*(-1)^(3-posinlink);


             normplane= links(linkid,6:8);

            nc(n,:)=nc(n,:)+ normplane;
            
            if l(linkid)>0
                fsegn0 =fseg(linkid,3*(posinlink-1)+(1:3));
                fn(n,:) = fn(n,:) + fsegn0; 
            end
            
        end
        nc(n,:)=nc(n,:)/norm(nc(n,:));
        fnc(n,:)=dot(fn(n,:),nc(n,:));
        fn_c=fnc.*nc;
        linvec(n,:)=linvec(n,:)/norm(linvec(n,:));
        vt(n,:)=linel(n)*linvec(n,:);


    end
    
    %%%assemble the matrix
    m=1; % the index of degrees
    z=zeros(nseg+nnode,1); % track the index of the degrees;
    loc_n0=0;
    loc_n1=0;
    loc_j=0;
    k_e=zeros(3,3);
    
    for i=1:nseg
        n0=links(i,1);
        n1=links(i,2);
       
        rt=rn(n1,1:3)-rn(n0,1:3);
        l(i)=norm(rt);
        burg=links(i,3:5);
        
        normplane= links(i,6:8);
        
        an0= nc(n0,:)*normplane';
        an1=nc(n1,:)*normplane';
        
        if l(i)>0
            
            theta=dot(rt,burg)/norm(burg)/l(i);
            theta=acos(theta);
            sintheta=abs(sin(theta));
        else
            sintheta=0;
        end
        
       %%% find the location of n0 n1 and middle node in the matrix; 
        if ismember(n0,z)
            loc_n0=find(z==n0);
        else 
            z(m)=n0;
            loc_n0=m;
            m=m+1;
        end
        z(m)=-i;
        loc_j=m;
        m=m+1;
        
        if ismember(n1,z)
            loc_n1=find(z==n1);
        else
            z(m)=n1;
            loc_n1=m;
            m=m+1;
        end

        k_e = l(i)*(bur*1e-9)/Dc* [23/960*(an0*l(i)*sintheta)^2*(bur*1e-9)^4,      17/960*(l(i)*sintheta)^2*an0*an1*(bur*1e-9)^4,   1/24*(an0*l(i)*sintheta)^1*(bur*1e-9)^2;
                                   17/960*(l(i)*sintheta)^2*an0*an1*(bur*1e-9)^4,  23/960*(l(i)*sintheta*an1)^2*(bur*1e-9)^4,       -1/24*(an1*l(i)*sintheta)^1*(bur*1e-9)^2;
                                   1/24*(an0*l(i)*sintheta)^1*(bur*1e-9)^2,        -1/24*(an1*l(i)*sintheta)^1*(bur*1e-9)^2,         1                           ]; 
        B_e = -[3/8*l(i)*sintheta*an0*(bur*1e-9)^2, 1/8*l(i)*sintheta*an1*(bur*1e-9)^2, 1;
               1/8*l(i)*sintheta*an0*(bur*1e-9)^2, 3/8*l(i)*sintheta*an1*(bur*1e-9)^2, -1];

                     
       if sintheta < 0.0872  %  cutt off angle 5;
           k_e=k_e*1e8;
           B_e=B_e*1e8;
        end
         
        k(loc_n0, loc_n0) = k(loc_n0, loc_n0) + k_e(1,1);
        k(loc_n0, loc_n1) = k(loc_n0, loc_n1) + k_e(1,2);
        k(loc_n0, loc_j)  = k(loc_n0, loc_j)  + k_e(1,3);
        k(loc_n1, loc_n0) = k(loc_n1, loc_n0) + k_e(2,1);
        k(loc_n1, loc_n1) = k(loc_n1, loc_n1) + k_e(2,2);
        k(loc_n1, loc_j)  = k(loc_n1, loc_j)  + k_e(2,3);
        k(loc_j,  loc_n0) = k(loc_j, loc_n0)  + k_e(3,1);
        k(loc_j,  loc_n1) = k(loc_j, loc_n1)  + k_e(3,2);
        k(loc_j,  loc_j)  = k(loc_j, loc_j)   + k_e(3,3);    
        
        C(n0,loc_n0) =  C(n0,loc_n0) + B_e(1,1);
        C(n0,loc_n1) =  C(n0,loc_n1) + B_e(1,2);
        C(n0,loc_j)  =  C(n0,loc_j)  + B_e(1,3);
        C(n1,loc_n0) =  C(n1,loc_n0) + B_e(2,1);
        C(n1,loc_n1) =  C(n1,loc_n1) + B_e(2,2);
        C(n1,loc_j)  =  C(n1,loc_j)  + B_e(2,3);
        
        Fnc(loc_n0) = fnc(n0);
        Fnc(loc_n1) = fnc(n1);
        Fnc(loc_j)  = 0;
        
    end
 
    
    
    
    % reshape the matrix
    lamla=zeros(nnode,nnode);
    K=[k C';C lamla];
    F=[Fnc*(bur*1e-9);zeros(nnode,1)];
    
    
    [vc,R]=linsolve(K,F); %% solve K*vc=F to get vc, R should be <<0 to ensure the accuracy of the solution;
 %    spy(K)  % show the structure of K
    

    %% rearrange the order of the degrees 

    for n=1:nnode
%        nnn=nc(n,:)/norm(nc(n,:));

       v_c(n,:)= vc(find(z==n))*nc(n,:);
       
     
    end

    
    rnmax=size(rn,2);
    for i=1:size(rn,1)
        if rn(i,rnmax)~=0
            v_c(i,1:3)=0;
        end
    end
    
    
end





