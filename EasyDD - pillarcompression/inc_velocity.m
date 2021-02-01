function [vn_g]=inc_velocity(vn_g,rn,links,connectivity,nodelist,conlist)

global RR QQ

angle=zeros(size(RR,2),1);

L0 = size(nodelist,1);

if L0==0
    %----calaulate vn_g for all nodes
    
    for i=1:size(links,1)
        
        P1=rn(links(i,1),1:3);
        P2=rn(links(i,2),1:3);
        
        for n=1:size(RR,2)
            r1=P1-P2;
            r2=QQ(n,1:3)-P2;
            r3=QQ(n,1:3)-P1;
            
            angle(n)=dot(r1,r2)*dot(-r1,r3);
            
            % plot3(rn([P1,P2],1),rn([P1,P2],2),rn([P1,P2],3));
            
            if angle(n)>=0  %acute trangle
                
                distance=norm(cross((P1-P2), (P2-QQ(n,:))))/norm(P1-P2);
                
                if distance <= RR(n)
                    
                    if dot(r3,vn_g(links(i,1),1:3))>=0
                        vn_g(links(i,1),1:3)=0;
                    end
                    
                    if dot(r2,vn_g(links(i,2),1:3))>=0
                        vn_g(links(i,2),1:3)=0;
                    end
                end
                
            else
                if norm(r3) <= RR(n) && dot(r3,vn_g(links(i,1),1:3))>=0
                    vn_g(links(i,1),1:3)=0;
                end
                
                if norm(r2) <= RR(n)&& dot(r2,vn_g(links(i,2),1:3))>=0
                    vn_g(links(i,2),1:3)=0;
                end
                
            end
            
        end
    end
    
else
    %----only calculate vn_g of the nodelist
    for j=1:L0
        n1=nodelist(j);
        numNbrs=conlist(j,1);
        for i=1:numNbrs
            nn=conlist(j,i+1);
            linkid=connectivity(n1,2*nn);
            posinlink=connectivity(n1,2*nn+1);
            n2=links(linkid,3-posinlink);
            
            P1=rn(n1,1:3);
            P2=rn(n2,1:3);
            
            for n=1:size(RR,2)
                r1=P1-P2;
                r2=QQ(n,1:3)-P2;
                r3=QQ(n,1:3)-P1;
                
                angle(n)=dot(r1,r2)*dot(-r1,r3);
                
                % plot3(rn([P1,P2],1),rn([P1,P2],2),rn([P1,P2],3));
                
                if angle(n)>=0  %acute trangle
                    
                    distance=norm(cross((P1-P2), (P2-QQ(n,:))))/norm(P1-P2);
                    
                    if distance <= RR(n)
                        
                        if dot(r3,vn_g(j,1:3))>=0
                            vn_g(j,1:3)=0;
                        end
                        
                        if dot(r2,vn_g(j,1:3))>=0
                            vn_g(j,1:3)=0;
                        end
                    end
                    
                else
                    if norm(r3) <= RR(n) && dot(r3,vn_g(j,1:3))>=0
                        vn_g(j,1:3)=0;
                    end
                    
                    if norm(r2) <= RR(n)&& dot(r2,vn_g(j,1:3))>=0
                        vn_g(j,1:3)=0;
                    end
                    
                end
                
            end
            
        end
    end
end


