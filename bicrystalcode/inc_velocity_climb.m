function [nearinc_flag]=inc_velocity_climb(vn_c,rn,links,nearinc_flag)

global RR QQ
angle=zeros(size(RR,2),1);

for i=1:size(links,1)
    
    
    P1=rn(links(i,1),1:3);
    P2=rn(links(i,2),1:3);
    
    for n=1:size(RR,2)
        r1=P1-P2;
        r2=QQ(n,1:3)-P2;
        r3=QQ(n,1:3)-P1;
        
        angle(n)=dot(r1,r2)*dot(-r1,r3);
        
        %      plot3(rn([P1,P2],1),rn([P1,P2],2),rn([P1,P2],3));
        
        if angle(n)>=0  %acute trangle
            
            distance=norm(cross((P1-P2), (P2-QQ(n,:))))/norm(P1-P2);
            
            if distance <= RR(n)
                
                if dot(r3,vn_c(links(i,1),1:3))>0
                    nearinc_flag(links(i,1))=1;
                end
                
                if dot(r2,vn_c(links(i,2),1:3))>0
                    nearinc_flag(links(i,2))=1;
                end
            end
            
        else
            if norm(r3) <= RR(n) && dot(r3,vn_c(links(i,1),1:3))>=0
                nearinc_flag(links(i,1))=1;
            end
            
            if norm(r2) <= RR(n)&& dot(r2,vn_c(links(i,2),1:3))>=0
                 nearinc_flag(links(i,2))=1;
            end
            
        end
        
    end

    
end
