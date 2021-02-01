function [vn_g]=inc_velocity_accelerate(vn_g,rn,links,connectivity,nodelist,conlist)


global RR QQ
% calculate the velocity near a inclusion
ninc=size(RR,2);
angle=zeros(ninc,1);
distance=angle;
flag = angle;

L0 = size(nodelist,1);

if L0==0
    %----calaulate vn_g for all nodes
    for i=1:size(links,1)
        
        n1=links(i,1);
        n2=links(i,2);
        P1=rn(n1,1:3);
        P2=rn(n2,1:3);
        
        r1=P1-P2;
        r2=P1-QQ;
        r3=P2-QQ;
        angle = (r2*r1').*(-r3*r1');
        distance = vecnorm(cross(r3, r2, 2)')/norm(r1);
        flag = (angle>=0)&(distance <= RR)';
        
        if any(flag==1)
            flag_ones=find(flag==1);
            inward1 = r2(flag_ones,:)*vn_g(n1,1:3)';
            inward2 = r3(flag_ones,:)*vn_g(n2,1:3)';
            vn_g(n1,1:3)= (1-any(inward1<=0))*vn_g(n1,1:3);
            vn_g(n2,1:3)= (1-any(inward2<=0))*vn_g(n2,1:3);
            
            
        else
            continue
        end
        
    end
    
else
    %----only calculate vn_g of the nodelist
    for j=1:L0
        n1=nodelist(j);
        numNbrs=conlist(j,1);
        for n=1:numNbrs
            nn=conlist(j,n+1);
            linkid=connectivity(n1,2*nn);
            posinlink=connectivity(n1,2*nn+1);
            n2=links(linkid,3-posinlink);
            
            P1=rn(n1,1:3);
            P2=rn(n2,1:3);
            
            r1=P1-P2;
            r2=P1-QQ;
            r3=P2-QQ;
            angle = (r2*r1').*(-r3*r1');
            distance = vecnorm(cross(r3, r2, 2)')/norm(r1);
            flag = (angle>=0)&(distance <= RR)';
            
            if any(flag==1)
                flag_ones=find(flag==1);
                inward1 = r2(flag_ones,:)*vn_g(j,1:3)';
                vn_g(j,1:3)= (1-any(inward1<=0))*vn_g(j,1:3);
           
            else
                continue
            end
        end
    end
    
    
end
end
