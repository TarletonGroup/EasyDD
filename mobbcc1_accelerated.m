function [vn,fn] = mobbcc1_accelerated(fseg,rn,links,connectivity,nodelist,conlist)
%mobility law function (model: BCC0)
global Bscrew Bedge Beclimb Bline

%numerical tolerance
eps=1e-12;

% length of the nodelist for which the velocity will be calculated
L1=size(nodelist,1);
% if no nodelist is given then the nodelist becomes the whole node population
% this portion of the code sets up that nodelist along with the connlist
% that contains all of the nodal connections
if L1==0
    L1=size(rn,1);
    nodelist=linspace(1,L1,L1)';
    [L2,L3]=size(connectivity);
    conlist=zeros(L2,(L3-1)/2+1);
    conlist(:,1)=connectivity(:,1);
    for i=1:L2
        connumb=conlist(i,1);
        conlist(i,2:connumb+1)=linspace(1,connumb,connumb);
    end
end

[fx,fy,fz,B11,B22,B33,B12,B13,B23]=mobbcc1mex(fseg,rn,links,connectivity,nodelist,conlist,Beclimb,Bline,Bscrew,Bedge);
fn=[fx,fy,fz];

for n=1:L1
    Btotal=[B11(n), B12(n), B13(n); B12(n), B22(n), B23(n); B13(n), B23(n), B33(n)];
    
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
    vn(n,:)=(evec*eval*evec'*fvec)';  % calculate the velocity 
else
    vn(n,:)=(Btotal\fn(n,:)')';                 % Btotal was wellconditioned so just take the inverse
end

    
%    if numNbrs==2
%        ii=conlist(n,2);                                                                      
%        n1=links(connectivity(n0,2*ii),3-connectivity(n0,2*ii+1));
%        ii=conlist(n,3);                                                                      
%        n2=links(connectivity(n0,2*ii),3-connectivity(n0,2*ii+1));
%        rt=rn(n1,1:3)-rn(n2,1:3);
%        L=norm(rt);
%        linedir=rt./L;
%        vn(n,:)=((eye(3)-linedir'*linedir)*vn(n,:)')';
%    end
end