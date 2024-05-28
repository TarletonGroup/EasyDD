function [vn_g]=touchboundary_velocity(vn_g,rn,links,connectivity,nodelist,conlist)

global SimBox rntol rmax

dim=1;%in x direction
loc_boundary=SimBox(1,dim)/2;% location of the boundary is x=dx/2

L0 = size(nodelist,1)
% norm_left=zeros(1,3);
nodelist

if L0==0
    for i=1:size(rn,1)
        
        if abs(rn(i,dim)-loc_boundary)<= rntol
            
            
            if rn(i,dim)< loc_boundary % on the left side of the boundary
                %             norm_left=norm_left(1,dim)-1;
                if dot(vn_g(i,dim),-1) < 0
                    vn_g(i,dim)=0;
                end
                
            elseif rn(i,dim)> loc_boundary
                if dot(vn_g(i,dim),1) < 0
                    vn_g(i,dim)=0;
                end
                
             elseif abs(rn(i,dim)-loc_boundary) == eps   
                 vn_g(i,1:3)=0;
            end
        end
        
    end

else
    for i=1:L0
        n1=nodelist(i);
        if abs(rn(n1,dim)-loc_boundary)<= rntol
            
            if rn(n1,dim)< loc_boundary % on the left side of the boundary
                %             norm_left=norm_left(1,dim)-1;
                if dot(vn_g(i,dim),-1) < 0
                    vn_g(i,dim)=0;
                end
                
            elseif rn(n1,dim)> loc_boundary
                if dot(vn_g(i,dim),1) < 0
                    vn_g(i,dim)=0;
                end
                
             elseif abs(rn(n1,dim)-loc_boundary) == eps   
                 vn_g(i,1:3)=0;
            end
        end        
        
    end
end

% if any(rn(:,1)> loc_boundary)
%     disp('nodes move accross the boundary! See touchboundary_velocity.m');
%     pause
% end





