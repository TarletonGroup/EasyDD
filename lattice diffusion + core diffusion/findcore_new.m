
function [gammaMiu_num,gammaMiu_el,gammaMiu_area] = findcore_new(rn,links,xnodes,wx,wy,wz,mx,my,mz,mel,nc)


% global bur
nseg=size(links,1);
planenorm=zeros(nseg,3);
burg=planenorm;
eps=1e-3;

% center of each element
xc=zeros(mel,3);
for n=1:mel
    for mm=1:3
        xc(n,mm)= 0.5*(xnodes(nc(n,1),mm) + xnodes(nc(n,7),mm)); %xc is center of element p
    end
end

%calculate the element number where the dislocation node located
gammaMiu_num=zeros(nseg,1);
gammaMiu_el=zeros(nseg,6);
gammaMiu_area=zeros(nseg,12);
for i=1:nseg
    n0=links(i,1);
    n1=links(i,2);
    
%     rt=-rn(n1,1:3)+ rn(n0,1:3);
    
    i0=max(ceil(rn(n0,1)/wx),1);
    j0=max(ceil(rn(n0,2)/wy),1);
    k0=max(ceil(rn(n0,3)/wz),1);
    
    p0=i0+(k0-1)*mx+(j0-1)*mx*mz;
    
    i1=max(ceil(rn(n1,1)/wx),1);
    j1=max(ceil(rn(n1,2)/wy),1);
    k1=max(ceil(rn(n1,3)/wz),1);
    
    p1=i1+(k1-1)*mx+(j1-1)*mx*mz;
    
    if p0==p1 %in a same element
        gammaMiu_el(i,1)= p0;
        gammaMiu_num(i)=1; %number of element invovled with seg i
        gammaMiu_area(i,1:2)=[-1 1]; %s[-1,1]
    else
        intersect_node=intersect(nc(p1,:),nc(p0,:));
        if size(intersect_node,2)==4 %two elements next to each other
            
            gammaMiu_el(i,1:2)=[p0 p1]; %two element involved
            gammaMiu_num(i)=2;
            
            %calculate the intersection point
            vec1=xnodes(intersect_node(2),:)- xnodes(intersect_node(1),:);
            vec2=xnodes(intersect_node(3),:)- xnodes(intersect_node(2),:);
            plan_n=cross(vec1,vec2); %normal of the commmon plan
            
            vec3=xnodes(intersect_node(1),:)-rn(n0,1:3);
            vec4=xnodes(intersect_node(1),:)-rn(n1,1:3);
            l1=vec3*plan_n';
            l2=vec4*plan_n';
            s=norm(l1)/(norm(l1)+norm(l2));
            gammaMiu_area(i,1:2)=[-1 2*s-1]; %p0: s[-1,s-1]; p1:s[s-1,1]
            gammaMiu_area(i,3:4)=[2*s-1 1];
            
        elseif size(intersect_node,2)==2 %two elements share a common ede
            
            %check the direction of the coomon edge:along x, y or z;
            vec1=xnodes(intersect_node(1),:) - xnodes(intersect_node(2),:);
            [m p]=max(abs(vec1));
            xyz1=rn(n1,1:3)-xnodes(intersect_node(1),1:3);
            xyz0=rn(n0,1:3)-xnodes(intersect_node(1),1:3);
            if p==1 %vec1 is along x
                %project into y-z plan
                z1z0=abs(xyz1(3))/abs(xyz0(3));
                y1y0=abs(xyz1(2))/abs(xyz0(2));
                if z1z0 > y1y0
                    k2=max(ceil((rn(n0,3)+ sign(rn(n1,3)-rn(n0,3))*wz)/wz),1);
                    p2=i0+(k2-1)*mx+(j0-1)*mx*mz;
                    
                    gammaMiu_el(i,1:3)=[p0 p2 p1]; %3 element involved
                    gammaMiu_num(i)=3;
                    gammaMiu_area(i,1:2)=[-1  2/(1+z1z0)-1];
                    gammaMiu_area(i,3:4)=[2/(1+z1z0)-1 1-2*y1y0/(1+y1y0)];
                    gammaMiu_area(i,5:6)=[1-2*y1y0/(1+y1y0)  1];
                    
                elseif abs(z1z0-y1y0)<=eps
                    gammaMiu_el(i,1:2)=[p0 p1]; %2 element involved
                    gammaMiu_num(i)=2;
                    s=1/(1+z1z0);
                    gammaMiu_area(i,1:2)=[-1  2*s-1];
                    gammaMiu_area(i,3:4)=[2*s-1  1];
                else
                    k2=max(ceil((rn(n1,3)+ sign(rn(n0,3)-rn(n1,3))*wz)/wz),1);
                    p2=i1+(k2-1)*mx+(j1-1)*mx*mz;
                    
                    gammaMiu_el(i,1:3)=[p0 p2 p1]; %3 element involved
                    gammaMiu_num(i)=3;
                    gammaMiu_area(i,1:2)=[-1  2/(1+y1y0)-1];
                    gammaMiu_area(i,3:4)=[2/(1+y1y0)-1 2*z1z0/(1+z1z0)-1];
                    gammaMiu_area(i,5:6)=[2*z1z0/(1+z1z0)-1  1];
                end
                
            elseif p==2  %vec1 is along y
                %project into x-z plan
                z1z0=abs(xyz1(3))/abs(xyz0(3));
                x1x0=abs(xyz1(1))/abs(xyz0(1));
                if z1z0 > x1x0
                    k2=max(ceil((rn(n0,3)+ sign(rn(n1,3)-rn(n0,3))*wz)/wz),1);
                    p2=i0+(k2-1)*mx+(j0-1)*mx*mz;
                    
                    gammaMiu_el(i,1:3)=[p0 p2 p1]; %3 element involved
                    gammaMiu_num(i)=3;
                    gammaMiu_area(i,1:2)=[-1  2/(1+z1z0)-1];
                    gammaMiu_area(i,3:4)=[2/(1+z1z0)-1 1-2*x1x0/(1+x1x0)];
                    gammaMiu_area(i,5:6)=[1-2*x1x0/(1+x1x0)  1];
                    
                elseif abs(z1z0-x1x0)<=eps
                    gammaMiu_el(i,1:2)=[p0 p1]; %2 element involved
                    gammaMiu_num(i)=2;
                    s=1/(1+z1z0);
                    gammaMiu_area(i,1:2)=[-1  2*s-1];
                    gammaMiu_area(i,3:4)=[2*s-1  1];
                else
                    k2=max(ceil((rn(n1,3)+ sign(rn(n0,3)-rn(n1,3))*wz)/wz),1);
                    p2=i1+(k2-1)*mx+(j1-1)*mx*mz;
                    
                    gammaMiu_el(i,1:3)=[p0 p2 p1]; %3 element involved
                    gammaMiu_num(i)=3;
                    gammaMiu_area(i,1:2)=[-1  2/(1+x1x0)-1];
                    gammaMiu_area(i,3:4)=[2/(1+x1x0)-1 1-2*z1z0/(1+z1z0)];
                    gammaMiu_area(i,5:6)=[1-2*z1z0/(1+z1z0)  1];
                end
            elseif p==3 % along z
                %project into x-y plan
                y1y0=abs(xyz1(2))/abs(xyz0(2));
                x1x0=abs(xyz1(1))/abs(xyz0(1));
                if y1y0 > x1x0
                    j2=max(ceil((rn(n0,2)+ sign(rn(n1,2)-rn(n0,2))*wy)/wy),1);
                    p2=i0+(k0-1)*mx+(j2-1)*mx*mz;
                    
                    gammaMiu_el(i,1:3)=[p0 p2 p1]; %3 element involved
                    gammaMiu_num(i)=3;
                    gammaMiu_area(i,1:2)=[-1  2/(1+y1y0)-1];
                    gammaMiu_area(i,3:4)=[2/(1+y1y0)-1 1-2*x1x0/(1+x1x0)];
                    gammaMiu_area(i,5:6)=[1-2*x1x0/(1+x1x0)  1];
                    
                elseif abs(x1x0-y1y0)<=eps
                    gammaMiu_el(i,1:2)=[p0 p1]; %2 element involved
                    gammaMiu_num(i)=2;
                    s=1/(1+y1y0);
                    gammaMiu_area(i,1:2)=[-1  2*s-1];
                    gammaMiu_area(i,3:4)=[2*s-1  1];
                else
                    j2=max(ceil((rn(n1,2)+ sign(rn(n0,2)-rn(n1,2))*wy)/wy),1);
                    p2=i1+(k1-1)*mx+(j2-1)*mx*mz;
                    
                    gammaMiu_el(i,1:3)=[p0 p2 p1]; %3 element involved
                    gammaMiu_num(i)=3;
                    gammaMiu_area(i,1:2)=[-1  2/(1+x1x0)-1];
                    gammaMiu_area(i,3:4)=[2/(1+x1x0)-1 1-2*y1y0/(1+y1y0)];
                    gammaMiu_area(i,5:6)=[1-2*y1y0/(1+y1y0)  1];
                end
            end
            
%           elseif size(intersect_node,2)==1 %two elements share a common node,at least 2 elements are invovled
%             pause  
        else
            if size(intersect_node,2)==0 || size(intersect_node,2)==1 %two elements do not share any node,at leat 3 element are invovled
            % center of each element
            distance_c=sqrt(wx^2+wy^2+wz^2)/2;
            gammaMiu_num(i)=2;
            gammaMiu_el(i,1)=p0;
            
            xc=zeros(mel,3);
            for n=1:mel
                for mm=1:3
                    xc(n,mm)= 0.5*(xnodes(nc(n,1),mm) + xnodes(nc(n,7),mm)); %xc is center of element p
                end
                
            end
            r1=-xc(:,1:3)+ rn(n0,1:3);
            r2=-xc(:,1:3)+ rn(n1,1:3);
            r12=rn(n0,1:3)-rn(n1,1:3);
            l=norm(r12);
            angel12=(r1*r12' >=1e-6).*(-r2*r12' >=1e-6);
            inbetween=find(angel12==1);
            for k=1: size(inbetween,1)
                n=inbetween(k);
                distance=norm(cross(r1(n,:),r2(n,:)))/l;
                if distance <= distance_c && n ~= p0 && n~= p1
                    gammaMiu_el(i,gammaMiu_num(i))=n;
                    gammaMiu_num(i)= gammaMiu_num(i) +1;
                    
                    %calculate the integration area (s0 s1)
                    if norm(cross(r12,[0 0 1]))>= 1e-6 %r12 is not parellel to axis z, project n0-n1 on xy plane 5-6-2-1
                        nr12=[-r12(2) r12(1)]; %norm of r12 on xy plane
                        n0_5=xnodes(nc(n,5),1:2)-rn(n0,1:2);
                        n0_6=xnodes(nc(n,6),1:2)-rn(n0,1:2);
                        n0_1=xnodes(nc(n,1),1:2)-rn(n0,1:2);
                        n0_2=xnodes(nc(n,2),1:2)-rn(n0,1:2);
                        
                        dotnr12_n05=dot(n0_5,nr12);
                        dotnr12_n06=dot(n0_6,nr12);
                        dotnr12_n01=dot(n0_1,nr12);
                        dotnr12_n02=dot(n0_2,nr12);
                        
%                         if dotnr12_n05*dotnr12_n06 < 1e-6 && dotnr12_n05*dotnr12_n01 < 1e-6
%                             %intersect with 5-6 and 5-1
%                         elseif dotnr12_n05*dotnr12_n06 < 1e-6 && dotnr12_n06*dotnr12_n02 < 1e-6
%                             %intersect with 5-6 and 6-2
%                         elseif dotnr12_n05*dotnr12_n06 < 1e-6 && dotnr12_n01*dotnr12_n02 < 1e-6
%                             %intersect with 5-6 and 1-2
%                         elseif dotnr12_n05*dotnr12_n01 < 1e-6 && dotnr12_n02*dotnr12_n06 < 1e-6
%                             %intersect with 5-1 and 2-6
%                         elseif dotnr12_n02*dotnr12_n06 < 1e-6 && dotnr12_n02*dotnr12_n01 < 1e-6
%                             %intersect with 2-6 and 1-2
%                         elseif dotnr12_n05*dotnr12_n01 < 1e-6 && dotnr12_n02*dotnr12_n01 < 1e-6
%                             %intersect with 5-1 and 1-2
%                         end
                        
                        if  dotnr12_n05*dotnr12_n06 < 1e-6 %intersect with 5-6
                            if dotnr12_n05*dotnr12_n01 < 1e-6 % 5-6/5-1
                               n1_5=xnodes(nc(n,5),1:2)-rn(n1,1:2);
                               
                               Y0Y1= norm(n0_5(2))/norm(n1_5(2));
                               X0X1= norm(n0_5(1))/norm(n1_5(1));
                               
                               aa=min(Y0Y1,X0X1);
                               bb=max(Y0Y1,X0X1);
                               
                               s0=(aa-1)/(aa+1);
                               s1=(bb-1)/(bb+1);
                               gammaMiu_area(i,(2*gammaMiu_num(i)-3) : (2*gammaMiu_num(i)-2) )=[s0 s1];
                               
                            elseif dotnr12_n02*dotnr12_n01 < 1e-6  %5-6/1-2
                                n1_5=xnodes(nc(n,5),1:2)-rn(n1,1:2);
                                n1_1=xnodes(nc(n,1),1:2)-rn(n1,1:2);
                                
                                Y05Y15=abs(n0_5(2))/abs(n1_5(2));
                                Y01Y11=abs(n0_1(2))/abs(n1_1(2));
                                
                                aa=min(Y05Y15,Y01Y11);
                                bb=max(Y05Y15,Y01Y11);
                                
                                s0=(aa-1)/(aa+1);
                                s1=(bb-1)/(bb+1);
                                gammaMiu_area(i,(2*gammaMiu_num(i)-3) : (2*gammaMiu_num(i)-2) )=[s0 s1];         
                                
                            elseif dotnr12_n02*dotnr12_n06 < 1e-6  %5-6/2-6
                                n1_6=xnodes(nc(n,6),1:2)-rn(n1,1:2);
                                
                                Y0Y1= norm(n0_6(2))/norm(n1_6(2));
                                X0X1= norm(n0_6(1))/norm(n1_6(1));
                                aa=min(Y0Y1,X0X1);
                                bb=max(Y0Y1,X0X1);
                                
                                s0=(aa-1)/(aa+1);
                                s1=(bb-1)/(bb+1);
                                gammaMiu_area(i,(2*gammaMiu_num(i)-3) : (2*gammaMiu_num(i)-2) )=[s0 s1];
                                
                            end
                        else
                            if dotnr12_n05*dotnr12_n01 < 1e-6 %5-1
                                if dotnr12_n02*dotnr12_n01 < 1e-6 %5-1/1-2
                                    n1_1=xnodes(nc(n,1),1:2)-rn(n1,1:2);
                                    Y0Y1= norm(n0_1(2))/norm(n1_1(2));
                                    X0X1= norm(n0_1(1))/norm(n1_1(1));
                                    aa=min(Y0Y1,X0X1);
                                    bb=max(Y0Y1,X0X1);
                                    
                                    s0=(aa-1)/(aa+1);
                                    s1=(bb-1)/(bb+1);
                                    gammaMiu_area(i,(2*gammaMiu_num(i)-3) : (2*gammaMiu_num(i)-2) )=[s0 s1];              
                                   
                                elseif dotnr12_n06*dotnr12_n02 < 1e-6 %5-1/6-2
                                    n1_5=xnodes(nc(n,5),1:2)-rn(n1,1:2);
                                    n1_2=xnodes(nc(n,2),1:2)-rn(n1,1:2);
                                    
                                    Y05Y15=abs(n0_5(2))/abs(n1_5(2));
                                    Y02Y12=abs(n0_2(2))/abs(n1_2(2));
                                    
                                    aa=min(Y05Y15,Y02Y12);
                                    bb=max(Y05Y15,Y02Y12);
                                    
                                    s0=(aa-1)/(aa+1);
                                    s1=(bb-1)/(bb+1);
                                    gammaMiu_area(i,(2*gammaMiu_num(i)-3) : (2*gammaMiu_num(i)-2) )=[s0 s1];
                                         
                                end
                                
                            elseif dotnr12_n01*dotnr12_n02 < 1e-6  %1-2
                                if dotnr12_n02*dotnr12_n06 < 1e-6  %1-2/2-6
                                    n1_2=xnodes(nc(n,2),1:2)-rn(n1,1:2);
                                    Y0Y1= norm(n0_2(2))/norm(n1_2(2));
                                    X0X1= norm(n0_2(1))/norm(n1_2(1));
                                    aa=min(Y0Y1,X0X1);
                                    bb=max(Y0Y1,X0X1);
                                    
                                    s0=(aa-1)/(aa+1);
                                    s1=(bb-1)/(bb+1);
                                    gammaMiu_area(i,(2*gammaMiu_num(i)-3) : (2*gammaMiu_num(i)-2) )=[s0 s1];
                                    
                                end
                            end
                        end                    
                        
                    end    
                    
                end
            end
            
            gammaMiu_el(i,gammaMiu_num(i))=p1;
            s0=min(gammaMiu_area(i,3:(2*gammaMiu_num(i)-2)));
            s1=max(gammaMiu_area(i,3:(2*gammaMiu_num(i)-2)));
            gammaMiu_area(i,1:2)=[-1 s0];
            gammaMiu_area(i,(2*gammaMiu_num(i)-1):2*gammaMiu_num(i))=[s1 1];
            end
            
        end
    end
end
 
    
end

% nseg=size(links,1);
% % planenorm1=zeros(nseg,3);
% % Fseg1=zeros(nseg,3);
% Miu_core1=zeros(nseg,1);
% % l1=Miu_core1;
% 
% tic
% rt1=rn(links(:,1),1:3)-rn(links(:,2),1:3);
% l1=vecnorm(rt1')';
% planenorm1=links(:,3:5);
% Fseg1=fseg(:,1:3)+fseg(:,4:6);
% Miu_core1=Fseg1.*planenorm1;
% Miu_core1=Miu_core1./l1;
% Miu_core1=sum(Miu_core1,2);
% 
% toc

%     my=mx;
%     mz=mx;
% ii=ceil(rn(n0,1)/wx);
% ii = max(ii,1);
% jj=ceil(rn(n0,2)/wy);
% jj =max(jj,1);
% kk=ceil(rn(n0,3)/wz);
% kk=max(kk,1);
% p = ii + (kk-1)*mx + (jj-1)*mx*mz;
% 
% ii=ceil(rn(n1,1)/wx);
% ii = max(ii,1);
% jj=ceil(rn(n1,2)/wy);
% jj =max(jj,1);
% kk=ceil(rn(n1,3)/wz);
% kk=max(kk,1);
% p1 = ii + (kk-1)*mx + (jj-1)*mx*mz;    
% 
% midnode=rn(n0,:)+rn(n1,:);
% 
% ii=ceil(midnode(1)/wx);
% ii = max(ii,1);
% jj=ceil(midnode(2)/wy);
% jj =max(jj,1);
% kk=ceil(midnode(3)/wz);
% kk=max(kk,1);
% pmid = ii + (kk-1)*mx + (jj-1)*mx*mz;