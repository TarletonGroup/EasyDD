function fseg = imageforce(MU,NU,a,rn,links,connectivity,fseg,Surpos,Btype);

rnlength = length(rn(:,1));
% add the image force
for n0=1:rnlength
    numNbrs=connectivity(n0,1);
  if ((Btype(3)) & (numNbrs==1))
    if  abs(rn(n0,3)-Surpos(3))<2*a | abs(rn(n0,3)+Surpos(3))<2*a
   %if ( rn(n0,3)==Surpos(3) ) & Btype(3) & (numNbrs==1)
        linkid = connectivity(n0,2);
        posinlink=connectivity(n0,3);
        
        n1=links(linkid,3-posinlink);
        rt=rn(n1,1:3)-rn(n0,1:3);                                                               % calculate the length of the link and its tangent line direction
        Lni=norm(rt);
        direc=rt/Lni;
        normal = cross( links(linkid,3:5),direc);
        nbeta = cross(direc,normal);
        beta = asin(norm(normal));
        
        rf(1)=rn(n0,1)-100;
        
        %% deal with some invalid link
        if links(linkid,7)==0 
            continue
        end
        rf(2)=links(linkid,6)*100/links(linkid,7) + rn(n0,2);
        %rf(3)=Surpos(3);
        if(rn(n0,3)>0)
            rf(3)=Surpos(3);
        else
             rf(3)=-Surpos(3);
        end
        inter=rn(n0,1:3)-rf(1:3);
        inter=inter/norm(inter);
        % get the angle between intection line and dislocation line
        project=dot(direc,inter);
        if (project<0)
            inter=-inter;
        end
        normal2=cross(direc,inter);
        nalfa=cross(normal2,direc);
        alfa=asin(norm(normal2));
        alfa=pi/2-alfa;
        
        cof1=abs( (1-NU*(cos(beta))^2)*tan(alfa) );
        cof2=abs(NU*sin(2*beta));
        image= (Lni/2)*MU/( 2*pi*(1-NU)*a )*( cof1*nalfa+cof2*nbeta );
        %image
        fseg( linkid,3*( posinlink-1)+1:3*( posinlink-1)+3 ) = fseg( linkid,3*( posinlink-1)+1:3*( posinlink-1)+3 ) + image;
    end
  end
    
  if ((Btype(1)) & (numNbrs==1))
    if  abs(rn(n0,1)-Surpos(1))<2*a | abs(rn(n0,1)+Surpos(1))<2*a
   %if ( rn(n0,3)==Surpos(3) ) & Btype(3) & (numNbrs==1)
        linkid = connectivity(n0,2);
        posinlink=connectivity(n0,3);
        
        n1=links(linkid,3-posinlink);
        rt=rn(n1,1:3)-rn(n0,1:3);                                                               % calculate the length of the link and its tangent line direction
        Lni=norm(rt);
        direc=rt/Lni;
        normal = cross( links(linkid,3:5),direc);
        nbeta = cross(direc,normal);
        beta = asin(norm(normal));
        
        rf(3)=rn(n0,3)-100;
        
        %% deal with some invalid link
        if links(linkid,7)==0 
            continue
        end
        rf(2)=links(linkid,8)*100/links(linkid,7) + rn(n0,2);
        
        %%%-----------------------
        %rf(3)=0;    %this part is used when the normal of slip plane is [0 0 1]
        %rf(2)=100;
        %%%______________________________
        
        if(rn(n0,1)>0)
            rf(1)=Surpos(1);
        else
             rf(1)=-Surpos(1);
        end
        inter=rn(n0,1:3)-rf(1:3);
        inter=inter/norm(inter);
        % get the angle between intection line and dislocation line
        project=dot(direc,inter);
        if (project<0)
            inter=-inter;
        end
        normal2=cross(direc,inter);
        nalfa=cross(normal2,direc);
        alfa=asin(norm(normal2));
        alfa=pi/2-alfa;
        
        cof1=abs( (1-NU*(cos(beta))^2)*tan(alfa) );
        cof2=abs(NU*sin(2*beta));
        image=(Lni/2)*MU/( 2*pi*(1-NU)*a )*( cof1*nalfa+cof2*nbeta );
        %image
        fseg( linkid,3*( posinlink-1)+1:3*( posinlink-1)+3 ) = fseg( linkid,3*( posinlink-1)+1:3*( posinlink-1)+3 ) + image;
    end
  end
  
  if ((Btype(2)) & (numNbrs==1))
    if  abs(rn(n0,2)-Surpos(2))<2*a | abs(rn(n0,2)+Surpos(2))<2*a
   %if ( rn(n0,3)==Surpos(3) ) & Btype(3) & (numNbrs==1)
        linkid = connectivity(n0,2);
        posinlink=connectivity(n0,3);
        
        n1=links(linkid,3-posinlink);
        rt=rn(n1,1:3)-rn(n0,1:3);                                                               % calculate the length of the link and its tangent line direction
        Lni=norm(rt);
        direc=rt/Lni;
        normal = cross( links(linkid,3:5),direc);
        nbeta = cross(direc,normal);
        beta = asin(norm(normal));
        
        rf(1)=rn(n0,1)-100;
        
        %% deal with some invalid link
        if links(linkid,8)==0 
            continue
        end
        rf(3)=links(linkid,6)*100/links(linkid,8) + rn(n0,3);
        %rf(2)=Surpos(2);
        if(rn(n0,2)>0)
            rf(2)=Surpos(2);
        else
             rf(2)=-Surpos(2);
        end
        inter=rn(n0,1:3)-rf(1:3);
        inter=inter/norm(inter);
        % get the angle between intection line and dislocation line
        project=dot(direc,inter);
        if (project<0)
            inter=-inter;
        end
        normal2=cross(direc,inter);
        nalfa=cross(normal2,direc);
        alfa=asin(norm(normal2));
        alfa=pi/2-alfa;
        
        cof1=abs( (1-NU*(cos(beta))^2)*tan(alfa) );
        cof2=abs(NU*sin(2*beta));
        image=(Lni/2)*MU/( 2*pi*(1-NU)*a )*( cof1*nalfa+cof2*nbeta );
        %image
        fseg( linkid,3*( posinlink-1)+1:3*( posinlink-1)+3 ) = fseg( linkid,3*( posinlink-1)+1:3*( posinlink-1)+3 ) + image;
    end
  end
end