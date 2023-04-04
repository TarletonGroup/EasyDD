function [rn,links]=inicon()
slipsys =[1 1 1 -1  1 0;
          1 1 1 -1  0 1;
          1 1 1  0  -1 1;
          -1 1 1  1 1 0;
          -1 1 1 1 0 1;
          -1 1 1 0 -1 1;
           1 -1 1 1 1 0;
           1 -1 1 1 0 1;
           1 -1 1 0 1 1;
           1 1 -1 -1 1 0;
           1 1 -1 1 0 1;
           1 1 -1 0 1 1];
lmax = 150;
lmin = 50;
elsize  =   (lmax+lmin)/2;
Simbox=[200 200 600];
Source_length = 200;
Source_devia = 0.01*Source_length;
nsegpsys = 8;
nsys=12;
nseg = nsegpsys*nsys;
elnum=0;
nodenum=0;
%%generate the random segment length 
length(1:nseg) = Source_length + Source_devia*randn(1,nseg);
%%generate the random start point for each segment
for i=1:3
    sego(1:nseg,i) = -Simbox(i) + 2*Simbox(i)*rand(nseg,1);
end

for j=1:nseg
    found = false;
    while ~found
        %% generat random line orintation
        fai(1)=(rand-0.5)*2*pi;
        lin(1)= cos(fai(1));

        if fai(1)>0
            lin(2)=sqrt(1-lin(1)^2);
        else
            lin(2)=-sqrt(1-lin(1)^2);
        end
        fai(3)=rand*pi;
        lin(3)=cos(fai(3));
        lin(1:3) = lin(1:3)/norm(lin(1:3));
        
        sege(j,1:3)=sego(j,1:3)+length(j)*lin(1:3);
        if (abs(sege(j,1))<Simbox(1)) & (abs(sege(j,2))<Simbox(2)) & (abs(sege(j,3))<Simbox(3))
            found = true;
        end
    end
    startp = sego(j,1:3);
    endp = sege(j,1:3);
    %% generate burger vertor and slip plane for each segment j
    if (mod(j,nsegpsys)==0)
       numbur = j/nsegpsys;
     else
       numbur =(j + nsegpsys - mod(j,nsegpsys))/nsegpsys;
     end
    burg = slipsys(numbur,4:6)/norm(slipsys(numbur,4:6));
    normplane = slipsys(numbur,1:3);
    if length(j)<elsize
        npart = 1;
    else
        npart   =   floor(length(j)/elsize);
    end
    %-----------------------------------------------------------------        
            for i=nodenum+1:nodenum+npart+1
                segnum=i-nodenum;
                alpha=(segnum-1)/npart;
                q(1:3)=(1-alpha).*startp + alpha.*endp;
                rn(i,1:3)=q(1:3);
                rn(i,4)=0;
                if(segnum==1)|(segnum==npart+1)
                    rn(i,4)=7;
                end
            end
           
            for i=elnum+1:elnum+npart
                links(i,1:2) = [i+j-1 i+j]; 
                links(i,3:5) = burg(1:3);
                links(i,6:8) = normplane(1:3);
            end
            elnum=elnum+npart;
            nodenum=nodenum+npart+1;
     %---------------------------------------------------------------------
end

plim=800;
viewangle=[135 45];
h=figure(1);
PlotBox(Simbox);
plotnodes(rn,links,plim); view(viewangle); xlim([-plim plim]);
ylim([-plim plim]); zlim([-plim plim]);
    
