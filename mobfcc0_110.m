function [vn,fn] = mobfcc0_110(fseg,rn,links,connectivity,nodelist,conlist)
%mobility law function (model: FCC0)

%Drag (Mobility) parameters (should be specified by Input file)
global Beclimb Bedge Bscrew Bline;
Bclimb=Beclimb;

%numerical tolerance
eps=1e-12;

%combinations of 111;
combinations111 = [1 1 1 ; -1 1 1 ; 1 -1 1 ; 1 1 -1 ; -1 -1 1 ; -1 1 -1 ; 1 -1 -1 ; -1 -1 -1];
rotationMatrix = [cosd(45) -sind(45) 0 ; sind(45) cosd(45) 0 ; 0 0 1];

for i=1:size(combinations111,1)
combinations111(i,:) = rotationMatrix*combinations111(i,:)'; %rotate based on orientation of slip planes.
combinations111(i,:) = combinations111(i,:) ./ norm( combinations111(i,:) ); %normalise
end

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
% now cycle through all of the nodes for which the velocity must be calculated

for n=1:L1
    n0=nodelist(n);                 %n0 is the nodeid of the nth node in nodelist
    numNbrs=conlist(n,1);           %numNbrs is the number of connections for node n0 in conlist
    fn(n,:)=zeros(1,3);             % initialize the total force and the total drag matrix
    Btotal=zeros(3,3);
    for i=1:numNbrs
        ii=conlist(n,i+1);                                                                      % connectionid for this connection
        linkid=connectivity(n0,2*ii);
        posinlink=connectivity(n0,2*ii+1);
        n1=links(linkid,3-posinlink);
        rt=rn(n1,1:3)-rn(n0,1:3);                                                               % calculate the length of the link and its tangent line direction
        L=norm(rt);
        if L>0.0
            fsegn0=fseg(linkid,3*(posinlink-1)+[1:3]);
            fn(n,:)=fn(n,:)+fsegn0; % nodeid for the node that n0 is connected to
            burgv=links(connectivity(n0,2*ii),3:5); % burgers vector of the link                                                           
            linedir=rt./L;
            nplane=links(linkid,6:8);
            nmag=norm(nplane);
            if nmag<eps
                % the normal plane is not defined try to define the normal plane
                Btotal=Btotal+(2.0*L).*((Bclimb).*eye(3)+(Bline-Bclimb).*(linedir'*linedir));
            else
                nplane=nplane./nmag;
                if any(sum(abs(repmat(nplane,8,1)-combinations111),2)/3<0.1) %nplane = <111> type
                    cth2=(linedir*burgv')^2/(burgv*burgv');                                                 % calculate how close to screw the link is
                    mdir=cross(nplane,linedir);
                    Bglide=1 / sqrt( 1 / Bedge^2 + ( 1 / Bscrew^2 - 1 / Bedge^2 ) * cth2);
                    Btotal=Btotal+(2.0*L).*( ( Bglide ).* ( mdir' * mdir ) + ( Bclimb ) .* ( nplane' * nplane ) +( Bline ).*( linedir' * linedir ) );
                else %nplane ~= <111> type. This is a special glide plane that is not of 111 type signifying a junction dislocation
                    Btotal=Btotal+(2.0*L).*((Bclimb).*eye(3)+(Bline-Bclimb).*(linedir'*linedir));
                end
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