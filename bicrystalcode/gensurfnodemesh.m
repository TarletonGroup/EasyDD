function [xsn,snc]=gensurfnodemesh(mx,my,mz,dx,dy,dz)

sx=dx/mx;
sy=dy/my;
sz=dz/mz;

[tx,ty,tz]=meshgrid(0:sx:dx,0:1:0,0:sz:dz);
Sfront=[tx(:) ty(:) tz(:)];

[tx,ty,tz]=meshgrid(0:1:0,sy:sy:dy,0:sz:dz);
Sleft=[tx(:) ty(:) tz(:)];

[tx,ty,tz]=meshgrid(dx:1:dx,sy:sy:dy,0:sz:dz);
Sright=[tx(:) ty(:) tz(:)];

[tx,ty,tz]=meshgrid(sx:sx:(dx-sx),sy:sy:dy,0:1:0);
tx=tx';
ty=ty';
Sbot=[tx(:) ty(:) tz(:)];

[tx,ty,tz]=meshgrid(sx:sx:(dx-sx),sy:sy:dy,dz:1:dz);
tx=tx';
ty=ty';
Stop=[tx(:) ty(:) tz(:)];

[tx,ty,tz]=meshgrid(sx:sx:(dx-sx),dy:1:dy,sz:sz:(dz-sz));
Sback=[tx(:) ty(:) tz(:)];

xsn=[Sfront;Sleft;Sright;Sbot;Stop;Sback];

Cfront=zeros(mx*mz,4);
for i=1:mz
    for j=1:mx
        ge=j+(i-1)*(mx);
        Cfront(ge,1)=j+(i-1)*(mx+1);
        Cfront(ge,2)=j+(i-1)*(mx+1)+1;
        Cfront(ge,3)=j+(i)*(mx+1);
        Cfront(ge,4)=j+(i)*(mx+1)+1;
    end
end


Cleft=zeros(my*mz,4);
for i=1:mz
    for j=1:my
        ge=j+(i-1)*(my);
        if mod(ge,my)==1
            Cleft(ge,1)=1+(i-1)*(mx+1);
            Cleft(ge,2)=(mx+1)*(mz+1)+j+(i-1)*my;
            Cleft(ge,3)=1+(i)*(mx+1);
            Cleft(ge,4)=(mx+1)*(mz+1)+j+i*my;
        else
            Cleft(ge,1)=(mx+1)*(mz+1)+(j-1)+(i-1)*my;
            Cleft(ge,2)=(mx+1)*(mz+1)+j+(i-1)*my;
            Cleft(ge,3)=(mx+1)*(mz+1)+(j-1)+i*my;
            Cleft(ge,4)=(mx+1)*(mz+1)+j+i*my;
        end
    end
end

Cright=zeros(my*mz,4);
for i=1:mz
    for j=1:my
       ge=j+(i-1)*my;
       if mod(ge,my)==1
          Cright(ge,1)=i*(mx+1);
          Cright(ge,2)=(mx+1)*(mz+1)+my*(mz+1)+j+(i-1)*my;
          Cright(ge,3)=(i+1)*(mx+1);
          Cright(ge,4)=(mx+1)*(mz+1)+my*(mz+1)+j+i*my;
       else
          Cright(ge,1)=(mx+1)*(mz+1)+my*(mz+1)+(j-1)+(i-1)*my;
          Cright(ge,2)=(mx+1)*(mz+1)+my*(mz+1)+j+(i-1)*my;
          Cright(ge,3)=(mx+1)*(mz+1)+my*(mz+1)+(j-1)+i*my;
          Cright(ge,4)=(mx+1)*(mz+1)+my*(mz+1)+j+i*my; 
       end
    end
end

Cbot=zeros(mx*my,4);
for i=1:my
   for j=1:mx
       gr=(i-1);
       ge=j+gr*mx;
       if mod(gr,my)==0
          if mod(ge,mx)==1
             Cbot(ge,1)=1;
             Cbot(ge,2)=2;
             Cbot(ge,3)=(mx+1)*(mz+1)+1;
             Cbot(ge,4)=(mx+1)*(mz+1)+2*my*(mz+1)+1;
          elseif mod(ge,mx)==0
             Cbot(ge,1)=mx;
             Cbot(ge,2)=mx+1;
             Cbot(ge,3)=(mx+1)*(mz+1)+2*my*(mz+1)+mx-1;
             Cbot(ge,4)=(mx+1)*(mz+1)+my*(mz+1)+1;
          else
             Cbot(ge,1)=j+(i-1)*mx;
             Cbot(ge,2)=j+(i-1)*mx+1;
             Cbot(ge,3)=(mx+1)*(mz+1)+2*my*(mz+1)+j+(i-1)*mx-1;
             Cbot(ge,4)=(mx+1)*(mz+1)+2*my*(mz+1)+j+(i-1)*mx; 
          end
       else
          if mod(ge,mx)==1
             Cbot(ge,1)=(mx+1)*(mz+1)+(i-1);
             Cbot(ge,2)=(mx+1)*(mz+1)+2*my*(mz+1)+j+(i-2)*(mx-1);
             Cbot(ge,3)=(mx+1)*(mz+1)+i;
             Cbot(ge,4)=(mx+1)*(mz+1)+2*my*(mz+1)+j+(i-1)*(mx-1);
          elseif mod(ge,mx)==0
             Cbot(ge,1)=(mx+1)*(mz+1)+2*my*(mz+1)+j-1+(i-2)*(mx-1);
             Cbot(ge,2)=(mx+1)*(mz+1)+my*(mz+1)+(i-1);
             Cbot(ge,3)=(mx+1)*(mz+1)+2*my*(mz+1)+j-1+(i-1)*(mx-1);
             Cbot(ge,4)=(mx+1)*(mz+1)+my*(mz+1)+i;
          else
             Cbot(ge,1)=(mx+1)*(mz+1)+2*my*(mz+1)+(j-1)+(i-2)*(mx-1);
             Cbot(ge,2)=(mx+1)*(mz+1)+2*my*(mz+1)+j+(i-2)*(mx-1);
             Cbot(ge,3)=(mx+1)*(mz+1)+2*my*(mz+1)+(j-1)+(i-1)*(mx-1);
             Cbot(ge,4)=(mx+1)*(mz+1)+2*my*(mz+1)+j+(i-1)*(mx-1); 
          end 
       end
   end
end

Ctop=zeros(mx*my,4);
for i=1:my
   for j=1:mx
       gr=(i-1);
       ge=j+gr*mx;
       if mod(gr,my)==0
          if mod(ge,mx)==1
             Ctop(ge,1)=1+(mx+1)*mz;
             Ctop(ge,2)=2+(mx+1)*mz;
             Ctop(ge,3)=(mx+1)*(mz+1)+1+my*mz;
             Ctop(ge,4)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*my+1;
          elseif mod(ge,mx)==0
             Ctop(ge,1)=(mx+1)*(mz+1)-1;
             Ctop(ge,2)=(mx+1)*(mz+1);
             Ctop(ge,3)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*my+mx-1;
             Ctop(ge,4)=(mx+1)*(mz+1)+my*(mz+1)+1+my*mz;
          else
             Ctop(ge,1)=j+(i-1)*mx+(mx+1)*mz;
             Ctop(ge,2)=j+(i-1)*mx+1+(mx+1)*mz;
             Ctop(ge,3)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*my+j+(i-1)*mx-1;
             Ctop(ge,4)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*my+j+(i-1)*mx; 
          end
       else
          if mod(ge,mx)==1
             Ctop(ge,1)=(mx+1)*(mz+1)+my*mz+i-1;
             Ctop(ge,2)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*my+j+(i-2)*(mx-1);
             Ctop(ge,3)=(mx+1)*(mz+1)+my*mz+i;
             Ctop(ge,4)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*my+j+(i-1)*(mx-1);
          elseif mod(ge,mx)==0
             Ctop(ge,1)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*my+j-1+(i-2)*(mx-1);
             Ctop(ge,2)=(mx+1)*(mz+1)+my*(mz+1)+my*mz+i-1;
             Ctop(ge,3)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*my+j-1+(i-1)*(mx-1);
             Ctop(ge,4)=(mx+1)*(mz+1)+my*(mz+1)+my*mz+i;
          else
             Ctop(ge,1)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*my+(j-1)+(i-2)*(mx-1);
             Ctop(ge,2)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*my+j+(i-2)*(mx-1);
             Ctop(ge,3)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*my+(j-1)+(i-1)*(mx-1);
             Ctop(ge,4)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*my+j+(i-1)*(mx-1); 
          end 
       end
   end
end

Cback=zeros(mx*mz,4);
for i=1:mz
   for j=1:mx
      gr=(i-1);
      ge=j+gr*mx;
      if mod(gr,mz)==0
         if mod(ge,mx)==1
            Cback(ge,1)=(mx+1)*(mz+1)+my;
            Cback(ge,2)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*(my-1)+1;
            Cback(ge,3)=(mx+1)*(mz+1)+2*my;
            Cback(ge,4)=(mx+1)*(mz+1)+2*my*(mz+1)+2*(mx-1)*my+1;
         elseif mod(ge,mx)==0
            Cback(ge,1)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*my;
            Cback(ge,2)=(mx+1)*(mz+1)+my*(mz+1)+my;
            Cback(ge,3)=(mx+1)*(mz+1)+2*my*(mz+1)+2*(mx-1)*my+(mx-1);
            Cback(ge,4)=(mx+1)*(mz+1)+my*(mz+1)+2*my; 
         else
            Cback(ge,1)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*(my-1)+j-1;
            Cback(ge,2)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*(my-1)+j;
            Cback(ge,3)=(mx+1)*(mz+1)+2*my*(mz+1)+2*(mx-1)*my+j-1;
            Cback(ge,4)=(mx+1)*(mz+1)+2*my*(mz+1)+2*(mx-1)*my+j; 
         end
      elseif mod(gr,mz)==(mz-1)
         if mod(ge,mx)==1
            Cback(ge,1)=(mx+1)*(mz+1)+my*mz;
            Cback(ge,2)=(mx+1)*(mz+1)+2*my*(mz+1)+2*(mx-1)*my+(mx-1)*(mz-2)+1;
            Cback(ge,3)=(mx+1)*(mz+1)+my*(mz+1);
            Cback(ge,4)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*my+(mx-1)*(my-1)+1;
         elseif mod(ge,mx)==0
            Cback(ge,1)=(mx+1)*(mz+1)+2*my*(mz+1)+2*(mx-1)*my+(mx-1)*(mz-1);
            Cback(ge,2)=(mx+1)*(mz+1)+my*(mz+1)+my*mz;
            Cback(ge,3)=(mx+1)*(mz+1)+2*my*(mz+1)+2*(mx-1)*my;
            Cback(ge,4)=(mx+1)*(mz+1)+2*my*(mz+1); 
         else
            Cback(ge,1)=(mx+1)*(mz+1)+2*my*(mz+1)+2*(mx-1)*my+j-1+(i-2)*(mx-1);
            Cback(ge,2)=(mx+1)*(mz+1)+2*my*(mz+1)+2*(mx-1)*my+j+(i-2)*(mx-1);
            Cback(ge,3)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*my+j-1+(mx-1)*(my-1);
            Cback(ge,4)=(mx+1)*(mz+1)+2*my*(mz+1)+(mx-1)*my+j+(mx-1)*(my-1); 
         end 
      else
         if mod(ge,mx)==1
            Cback(ge,1)=(mx+1)*(mz+1)+i*my;
            Cback(ge,2)=(mx+1)*(mz+1)+2*my*(mz+1)+2*(mx-1)*my+j+(i-2)*(mx-1);
            Cback(ge,3)=(mx+1)*(mz+1)+(i+1)*my;
            Cback(ge,4)=(mx+1)*(mz+1)+2*my*(mz+1)+2*(mx-1)*my+j+(i-1)*(mx-1);
         elseif mod(ge,mx)==0
            Cback(ge,1)=(mx+1)*(mz+1)+2*my*(mz+1)+2*(mx-1)*my+j-1+(i-2)*(mx-1);
            Cback(ge,2)=(mx+1)*(mz+1)+my*(mz+1)+i*my;
            Cback(ge,3)=(mx+1)*(mz+1)+2*my*(mz+1)+2*(mx-1)*my+j-1+(i-1)*(mx-1);
            Cback(ge,4)=(mx+1)*(mz+1)+my*(mz+1)+(i+1)*my; 
         else
            Cback(ge,1)=(mx+1)*(mz+1)+2*my*(mz+1)+2*(mx-1)*my+j-1+(i-2)*(mx-1);
            Cback(ge,2)=(mx+1)*(mz+1)+2*my*(mz+1)+2*(mx-1)*my+j+(i-2)*(mx-1);
            Cback(ge,3)=(mx+1)*(mz+1)+2*my*(mz+1)+2*(mx-1)*my+j-1+(i-1)*(mx-1);
            Cback(ge,4)=(mx+1)*(mz+1)+2*my*(mz+1)+2*(mx-1)*my+j+(i-1)*(mx-1); 
         end 
      end
   end
end

snc=[Cfront;Cleft;Cright;Cbot;Ctop;Cback];

end