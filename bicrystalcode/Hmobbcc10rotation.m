function [vn,fn] = Hmobbcc10rotation(fseg,rn,links,connectivity,nodelist,conlist,concentration)
%mobility law function (model: FCC0)

% [rn,links] = cross_slip_FCCrotate(fseg,rn,links,connectivity,nodelist,conlist);

global Bscrew Bedge Beclimb Bline rotationBCC

rotationBCC = eye(3);

rn0 = rn;
links0 = links;

combinations110 = 1/sqrt(2) * [1 1 0 ; -1 -1 0 ; 1 -1 0 ; -1 1 0 ; ...
    1 0 1 ; -1 0 -1 ; 1 0 -1 ; -1 0 1 ; ...
    0 1 1 ; 0 -1 -1 ; 0 1 -1 ; 0 -1 1 ];

%HY20181010: rotate rn (line direction) to crystal system
rn(:,1:3) = rn(:,1:3)*rotationBCC';
%HY20181010: the fseg and links also need to be rotated
fseg(:,1:3) = fseg(:,1:3)*rotationBCC';
fseg(:,4:6) = fseg(:,4:6)*rotationBCC';
links(:,3:5) = links(:,3:5)*rotationBCC';
links(:,6:8) = links(:,6:8)*rotationBCC';

%numerical tolerance
eps=1e-6;

global vdotl vnorm vnormvec

% vdotl = zeros(size(rn,1),1);
% vnorm = zeros(size(rn,1),1);
% vnormvec = zeros(size(rn,1),3);

links(find(abs(links)<eps)) = 0;

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

highlighti = 0;
linklater = 0;

linedir_ave = zeros(L1,3);
neighbours = zeros(L1,2);
vdotl = zeros(L1,1);
vnorm = zeros(L1,1);
vnormvec = zeros(L1,3);

for n=1:L1
    n0=nodelist(n);                 %n0 is the nodeid of the nth node in nodelist
    %     if rn(n0,end)==7
    %         continue
    %     end
    numNbrs=conlist(n,1);           %numNbrs is the number of connections for node n0 in conlist
    fn(n,:)=zeros(1,3);             % initialize the total force and the total drag matrix
    Btotal=zeros(3,3);
    linetemp = zeros(2,3);
    nplanetemp = zeros(numNbrs,3);
    nplaneX = zeros(numNbrs,3);
    %     neighbours = zeros(L1,2);
    for i=1:numNbrs
        ii=conlist(n,i+1);                                                                      % connectionid for this connection
        linkid=connectivity(n0,2*ii);
        if linkid>size(fseg,1)
%             disp('linkid>size(fseg,1)')
            continue
        end
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
%             nplaneX(i,:) = nplane;
            nplanetemp(i,:) = nplane;
            
            if nmag<eps
                % the normal plane is not defined try to define the normal plane
                Btotal=Btotal+(0.5*L).*((Beclimb).*eye(3)+(Bline-Beclimb).*(linedir'*linedir));
%                 disp('the normal plane is not defined')
                %                 pause
            else
                nplane=nplane./nmag;
                nplaneX(i,:) = nplane;
                nplanetemp(i,:) = nplane;
                if linedir*nplane'>1E-2 && rn(n0,4)~=67
                    % disp('linedir*nplane>eps')
                    highlighti = highlighti+1;
                    highlight(highlighti) = n0;
                    % pause
                    % linedir = linedir-(linedir*nplane')*nplane;
                    % pause
                end
                if any(sum(abs(repmat(nplane,12,1)-combinations110),2)/3<0.1) %nplane = <111> type
                    cth2=(linedir*burgv')^2/(burgv*burgv');                                                 % calculate how close to screw the link is
                    mdir=cross(nplane,linedir);
                    Bglide=1 / sqrt( 1 / Bedge^2 + ( 1 / Bscrew^2 - 1 / Bedge^2 ) * cth2);
                    Btotal=Btotal+(0.5*L).*( ( Bglide ).* ( mdir' * mdir ) + ( Beclimb ) .* ( nplane' * nplane ) +( Bline ).*( linedir' * linedir ) );
                    if numNbrs == 2
                        linetemp(i,:) = linedir;
                        neighbours(n,i) = n1;
                    end
                else %nplane ~= <111> type. This is a special glide plane that is not of 111 type signifying a junction dislocation
                    Btotal=Btotal+(0.5*L).*((Beclimb).*eye(3)+(Bline-Beclimb).*(linedir'*linedir));
%                     disp('climbing')
                end
            end
        end
    end
    
    if numNbrs == 2
        linedir_ave(n,:) = (linetemp(1,:)-linetemp(2,:))/sqrt(2);
    end
    
    if numNbrs == 2 && norm(linetemp(1,:)+linetemp(2,:))<1E-6
        %         Btotal=2*((0.5*L).*( ( Bglide ).* ( mdir' * mdir ) + ( Beclimb ) .* ( nplane' * nplane ) +( Bline ).*( linedir' * linedir )));
        %         Btotal=2*((0.5*L).*( ( Bglide ).* ( mdir' * mdir ) + ( Beclimb ) .* ( nplane' * nplane ) +( Bglide ).*( linedir' * linedir )));
    end
    
    %     if (numNbrs==2 & norm(nplaneX(1,:)-nplaneX(2,:))>1E-4 &...
    %         norm(nplaneX(1,:)+nplaneX(2,:))>1E-4 & norm(nplaneX(1,:))>0.1 &...
    %         norm(nplaneX(2,:))>0.1)
    %         linklater = linklater+1;
    %         linedirX = cross(nplaneX(1,:),nplaneX(2,:));
    %         linedirX = linedirX/norm(linedirX);
    %         Btotal = Btotal+(0.5*L).*((Beclimb).*eye(3)+(Bline-Beclimb).*(linedirX'*linedirX));
    % %         Btotal=BtotalX+(0.5*L).*((Beclimb).*eye(3)+(Bline-Beclimb).*(linedirX'*linedirX));
    % %         Btotal=(0.5*L).*((Beclimb).*eye(3)+(Bline-Beclimb).*(linedir'*linedir));
    %     end
    %     rcond(Btotal)
    if rcond(Btotal)<eps
        
%         %         if rn(n0,end)~=7
%         %         disp('rcond(Btotal)<eps')
%         %         end
%         
%         [evec,eval]=eig(Btotal);                    % find eigenvalues and eigen vectors of drag matrix
%         evalmax=eval(1,1);
%         eval=eval./evalmax;
%         fvec=fn(n,:)'./evalmax;
%         for i=2:3                                   % invert drag matrix and keep zero eigen values as zero
%             if eval(i,i)>eps
%                 eval(i,i)=1/eval(i,i);
%             else
%                 eval(i,i)=0.0d0;
%             end
%         end
%         vn(n,:)=(evec*eval*evec'*fvec)';  % calculate the velocity
%         %         vn(n,:)./mdir
%         %         dot(linedir,vn(n,:))
        Btotal = Btotal + eye(3)*max(max(abs(Btotal)))*10^(-6);
        vn(n,:)=(Btotal\fn(n,:)')';
    else
        vn(n,:)=(Btotal\fn(n,:)')';                 % Btotal was wellconditioned so just take the inverse
        %         vn(n,:)
        %         mdir
        %         vn(n,:)./mdir
        %         dot(linedir,vn(n,:))
    end
    
    %     if n==330
    %         pause
    %     end
    
    %     for kk=1:numNbrs
    %         vn(n,:) = vn(n,:)-dot(vn(n,:),nplanetemp(kk,:))*nplanetemp(kk,:);
    %     end
    
    %     vn(n,:) = vn(n,:)-dot(vn(n,:),nplane)*nplane;
    
    %     [vn(n,:),iflag,it]=Fp_SysLin(Btotal,fn(n,:)');
    %     if iflag==0
    %         disp('inaccurate vn')
    %         pause
    %     end
    
    %     vn(n,:) = vn(n,:)-dot(vn(n,:),nplane)*nplane;
    
    %     vdotl(n,1) = dot(linedir_ave(n,:),vn(n,:));
    %     vnorm(n,:) = vn(n,:) - vdotl(n,1)*linedir_ave(n,:);
    
    if numNbrs==2 && exist('nplane','var')
        if (norm(nplaneX(1,:)-nplaneX(2,:))>1E-4 &...
                norm(nplaneX(1,:)+nplaneX(2,:))>1E-4 & norm(nplaneX(1,:))>0.1 &...
                norm(nplaneX(2,:))>0.1)
            %         linklater = linklater+1;
            linedirX = cross(nplaneX(1,:),nplaneX(2,:));
            linedirX = linedirX/norm(linedirX);
            vn(n,:) = dot(vn(n,:),linedirX).*linedirX;
        else
            vn(n,:) = vn(n,:)-dot(vn(n,:),nplane)*nplane;%HY:weighted average of all the nplanes
        end
        
        linedir_ave(n,:) = linedir_ave(n,:)/norm(linedir_ave(n,:));
        vdotlvec(n,:) = dot(linedir_ave(n,:),vn(n,:))*linedir_ave(n,:);
        vdotl(n,1) = norm(vdotlvec(n,:));
        vnormvec(n,:) = vn(n,:) - vdotl(n,1)*linedir_ave(n,:);
        vnorm(n,1) = norm(vnormvec(n,:));
    end
end

%HY20180727: added by HY for HJunction problem to solve the infinite
%iteration, together with a modification in int_trapzoidal. Should not be taken as a permanent solution though.
% [row, col] = find(isnan(vn));
% vn(row, col) = 0;

[row, col] = find(rn(nodelist,end)==7);
vn(row, 1:3) = 0;
fn(row, 1:3) = 0;
vdotl(row, 1) = 0;
[row, col] = find(rn(nodelist,end)==8);
vn(row, 1:3) = 0;
fn(row, 1:3) = 0;
vdotl(row, 1) = 0;

% for kk=1:L1
%     if numNbrs==2 && norm(linedir_ave(kk,:))>eps && norm(linetemp(1,:)+linetemp(2,:))<1E-6
%         linedir_ave(kk,:) = linedir_ave(kk,:)/norm(linedir_ave(kk,:));
%         vdotlvec(kk,:) = dot(linedir_ave(kk,:),vn(kk,:))*linedir_ave(kk,:);
%         vdotl(kk,1) = norm(vdotlvec(kk,:));
%         vnormvec(kk,:) = vn(kk,:) - vdotl(kk,1)*linedir_ave(kk,:);
%         vnorm(kk,1) = norm(vnormvec(kk,:));
%     end
% end

% vdotl = abs(sum(linedir_ave.*vn,2));
% [row, col] = find(rn(nodelist,end)==7);
% vn(row, 1:3) = 0;
% fn(row, 1:3) = 0;

%HY20181010: rotate vn back to global system
vn = vn*rotationBCC;
%HY20181018: fn should also be rotated back!!!!
fn = fn*rotationBCC;

% if highlighti>0
% plim = 1.1314e+03;
% figure(111);
% highlightnodes([rn0(:,1:3),rn0(:,4)],links0,plim,highlight);
% drawnow
% % disp('highlighted nodes')
% end

%HY20180727: added by HY for HJunction problem to solve the infinite
%iteration, together with a modification in int_trapzoidal. Should not be taken as a permanent solution though.
[row, col] = find(isnan(vn));
vn(row, col) = 0;