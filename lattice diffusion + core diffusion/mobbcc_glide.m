% modified by gaoyuan 09-09-18

function [vn_g,fn] = mobbcc_glide(fseg,rn,links,connectivity,nodelist,conlist)

global Bscrew Bedge Beclimb Bline

vmax  = 3.9e11;        % unit:b/s,the max velocity is 100 m/s
eps   = 1e-12;

L0 = size(nodelist,1);

if L0==0
    L1 = size(rn,1);    % when mobility law is call form remesh fucnction, L1~=L0;
    nodelist = linspace(1,L1,L1)';
    [L2,L3] = size(connectivity);
    conlist = zeros(L2,(L3-1)/2+1);
    conlist(:,1) = connectivity(:,1);
    for i=1:L2
        connumb = conlist(i,1);
        conlist(i,2:connumb+1) = linspace(1,connumb,connumb);
    end
else
    L1=L0;
end

vn_g = zeros(L1,3);
fn=zeros(L1,3);

for n=1:L1
    n0=nodelist(n);                 %n0 is the nodeid of the nth node in nodelist
    numNbrs=conlist(n,1);           %numNbrs is the number of connections for node n0 in conlist
    % initialize the total force and the total drag matrix
    Btotal=zeros(3,3);
    for i=1:numNbrs
        ii=conlist(n,i+1);                                                                      % connectionid for this connection
        linkid=connectivity(n0,2*ii);
        posinlink=connectivity(n0,2*ii+1);
        n1=links(linkid,3-posinlink);
        rt=rn(n1,1:3)-rn(n0,1:3);                                                          % calculate the length of the link and its tangent line direction
        L=norm(rt);
        %fprintf('ii=%i, linkid=%i, n0=%i, n1=%i, L=%f \n',ii,linkid,n0,n1,L);
        if L>0.0
            fsegn0=fseg(linkid,3*(posinlink-1)+[1:3]);
            fn(n,:)=fn(n,:)+fsegn0; % nodeid for the node that n0 is connected to
            burgv=links(connectivity(n0,2*ii),3:5); % burgers vector of the link
            linedir=rt./L;
            
            
            costh2=(linedir*burgv')^2/(burgv*burgv'); % (lhat.bhat)^2 = cos^2(theta)   % calculate how close to screw the link is
            sinth2 = 1-costh2;
            Btotal=Btotal+(0.5*L).*((Bscrew).*eye(3)+(Bline-Bscrew).*(linedir'*linedir));   % build the drag matrix assuming that the dislocation is screw type
            
            if sinth2 > eps % not pure screw segment
                % correct the drag matrix for dislocations that are not screw type
                ndir=cross(burgv,linedir)./sqrt((burgv*burgv')*sinth2); % ndir = bxl/norm(bxl)
                mdir=cross(ndir,linedir);
                
                Bglide=1 / sqrt( (1 / Bedge^2) * sinth2 + ( 1 / Bscrew^2 ) * costh2); % Eqn (112) from Arsenlis et al 2007 MSMSE 15 553
                Bclimb=sqrt( (Beclimb^2 ) * sinth2 + ( Bscrew^2 ) * costh2);
                Btotal=Btotal+(0.5*L).*(( Bglide - Bscrew ).* ( mdir' * mdir ) + ( Bclimb - Bscrew ) .* ( ndir' * ndir ) );
            end
        end
    end
    
    
    if rcond(Btotal)<eps
        [evec,eval]=eig(Btotal);                    % find eigenvalues and eigen vectors of drag matrix
        evalmax=eval(1,1);
        eval=eval./evalmax;
        fvec=fn(n,:)'./evalmax;
        for i=2:3                                   % invert drag matrix and keep zero eigen values as zero
            if eval(i,i)>eps
                eval(i,i)=1/eval(i,i);
            else
                eval(i,i)=0.0d0;
            end
        end
        vn_g(n,:)=(evec*eval*evec'*fvec)';  % calculate the velocity
    else
        vn_g(n,:)=(Btotal\fn(n,:)')';                 % Btotal was wellconditioned so just take the inverse
    end
    
end

rnmax=size(rn,2);
for i=1:size(rn,1)
    if rn(i,rnmax)==4
        vn_g(i,1:3)=0;
    end
end

end





