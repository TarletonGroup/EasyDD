% modified by gaoyuan 09-09-18

function [vn,links,fn_g] = mobbcc_linear_full_glide_test(fseg,rn,links,connectivity,nodelist,conlist)

% vn  total velocity
% vn_c climb velocity
% vn_g glide velocity
% fn1 total nodal force

%mobility law function (model: FCC1)
%FCC1: velocity linear to force but remain orthogonal to glide plane
% if nodelist and conlist are empty then vn is the length of rn otherwise
% vn is the length of nodelist
global bur tem do_inclusion

vmax  = 3.9e11;        % unit:b/s,the max velocity is 100 m/s
tol=1e-7; %numerical tolerance
Bedge=1e0;

% drag coefficient for glide
Bscrew=1e0;
Bedge=1e0;
Bclimb=1e9;
Bline=1.0e-4*min(Bscrew,Bedge);

nodelist1=nodelist;
conlist1 =conlist;

%------- length of the nodelist for which the velocity will be calculated
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

% now cycle through all of the nodes for which the velocity must be calculated
vn_g=zeros(L1,3);  % nodal velocity
fn_g=zeros(L1,3); % nodal force

vn1=vn_g;   

rnmax=size(rn,2);

for n=1:L1
    n0=nodelist(n);
    numNbrs=conlist(n,1);
    nv=zeros(numNbrs,3);

    L=0;
    
    nplanetemp = zeros(numNbrs,3);
    nplaneX = zeros(numNbrs,3);
    for i=1:numNbrs
        ii=conlist(n,i+1);
        linkid=connectivity(n0,2*ii);
        posinlink=connectivity(n0,2*ii+1);
        n1=links(linkid,3-posinlink);
        rt=rn(n1,1:3)-rn(n0,1:3);
        Lni=norm(rt);
        L = L + Lni;
        if Lni>0
            fsegn0 =fseg(linkid,3*(posinlink-1)+[1:3]);
            fn_g(n,:)=fn_g(n,:)+fsegn0;
            nv(i,:)=links(linkid,6:8);
            linevec=rt/Lni;
            %             nplaneX(i,:) = nv(i,:);
            %             nplanetemp(i,:) = nplane;
        end
    end
    
    vn1(n,:)=fn_g(n,:);
    L=L/2;
    
    if (L==0)
        continue;
    end
    
    vn1(n,:)=vn1(n,:)/(L*Bedge); %(active in mobfcc1)
    
    if numNbrs==2
        if (norm(nv(1,:)-nv(2,:))>1E-4 &...
                norm(nv(1,:)+nv(2,:))>1E-4 & norm(nv(1,:))>0.1 &...
                norm(nv(2,:))>0.1)
            %         linklater = linklater+1;
            linedirX = cross(nv(1,:),nv(2,:));
            linedirX = linedirX/norm(linedirX);
            vn1(n,:) = dot(vn1(n,:),linedirX).*linedirX;
            %         else
            %             vn_g(n,:) = vn_g(n,:)-dot(vn_g(n,:),nplane)*nplane;%HY:weighted average of all the nplanes
        end
        %
        %         linedir_ave(n,:) = linedir_ave(n,:)/norm(linedir_ave(n,:));
        %         vdotlvec(n,:) = dot(linedir_ave(n,:),vn_g(n,:))*linedir_ave(n,:);
        %         vdotl(n,1) = norm(vdotlvec(n,:));
        %         vn_gormvec(n,:) = vn_g(n,:) - vdotl(n,1)*linedir_ave(n,:);
        %         vn_gorm(n,1) = norm(vn_gormvec(n,:));
    end
    
    
    %  orthogonalize vn to all gilde plane normal vectors
    
    %     if rn(n0,rnmax)~=0
    %            normal = nv(1,:);
    %         if rn(n0,rnmax)==-11
    %                n2n3=normal(2)^2+normal(3)^2;
    %                if n2n3~=0
    %                    vn1(n,1:3) =(vn1(n,3)*normal(2)-vn1(n,2)*normal(3))/n2n3*[0 -normal(3) normal(2)];
    %                else
    %                    vn1(n,1:3) = [0 0 0];
    %                end
    %         elseif rn(n0,rnmax)==-12
    %               n1n3=normal(1)^2+normal(3)^2;
    %                if n1n3~=0
    %                    vn1(n,1:3) =(vn1(n,3)*normal(1)-vn1(n,1)*normal(3))/n1n3*[-normal(3) 0 normal(1)];
    %                else
    %                    vn1(n,1:3) = [0 0 0];
    %                end
    %         elseif rn(n0,rnmax)==-13
    %                n1n2=normal(1)^2+normal(2)^2;
    %                if n1n2~=0
    %                    vn1(n,1:3) =(vn1(n,2)*normal(1)-vn1(n,1)*normal(2))/n1n2*[-normal(2) normal(1) 0 ];
    %                else
    %                    vn1(n,1:3) = [0 0 0];
    %                end
    % %         elseif rn(n0,rnmax)==4
    % %               vn1(n0,1:3) = vn1(n0,1:3)-vn1(n0,1:3)*;
    %         end
    %         if norm(vn1(n,:))>vmax
    %             vn1(n,:)=vmax*vn1(n,:)/norm(vn1(n,:));
    %         end
    % %         continue
    %     end
    
    for i=1:numNbrs
        for j=1:i-1,
            nv(i,:)=nv(i,:)-dot(nv(i,:),nv(j,:))*nv(j,:);
        end
        if(norm(nv(i,:))>eps)
            nv(i,:)=nv(i,:)/norm(nv(i,:));
        else
            nv(i,:)=0;
        end
    end
    for i=1:numNbrs
        vn1(n,:)=vn1(n,:)-dot(vn1(n,:),nv(i,:))*nv(i,:);
    end
    
    % ensure the velocity don't exceed 100m/s,because the Bedge increase
    % with velocity
    
    if norm(vn1(n,:))>vmax
        vn1(n,:)=vmax*vn1(n,:)/norm(vn1(n,:));
    end
    
    if rn(n,rnmax)==7
        %             vn1(n,1:3) = dot(vn1(n,1:3) ,links(linkid,3:5))*links(linkid,3:5);
        vn1(n,1:3) = [0 0 0];
        
    end
    
    
end

if do_inclusion ==1
 %   [vn]=inc_velocity_accelerate(vn1,rn,links,connectivity,nodelist1,conlist1 ); % calculate the dislocation glide mobility near inclusions
    [vn]=inc_velocity(vn1,rn,links,connectivity,nodelist1,conlist1 );
else
    vn=vn1;
end










