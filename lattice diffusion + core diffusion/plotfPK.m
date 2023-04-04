function plotfPK(rn,links,connectivity,plim,fseg)
%HY20180603: created by HY to plot PK forces on the dislocation line

L1=size(rn,1);
nodelist=linspace(1,L1,L1)';
[L2,L3]=size(connectivity);
conlist=zeros(L2,(L3-1)/2+1);
conlist(:,1)=connectivity(:,1);
for i=1:L2
    connumb=conlist(i,1);
    conlist(i,2:connumb+1)=linspace(1,connumb,connumb);
end
fn = zeros(L1,3);     
Hfn = zeros(L1,3); 
seglen = zeros(L1,1);
flag = zeros(L1,1);
for n=1:L1
    n0=nodelist(n);                 %n0 is the nodeid of the nth node in nodelist
    numNbrs=conlist(n,1);           %numNbrs is the number of connections for node n0 in conlist
    for i=1:numNbrs
        ii=conlist(n,i+1);                                                                      % connectionid for this connection
        linkid=connectivity(n0,2*ii);
        posinlink=connectivity(n0,2*ii+1);
        n1=links(linkid,3-posinlink);
        rt=rn(n1,1:3)-rn(n0,1:3);                                                          % calculate the length of the link and its tangent line direction
        L=norm(rt);
        seglen(n,1) = L;
        if L>0.0
            fsegn0=fseg(linkid,3*(posinlink-1)+[1:3]);
            fn(n,:)=fn(n,:)+fsegn0;  
%             Hfsegn0=Hfseg(linkid,3*(posinlink-1)+[1:3]);
%             Hfn(n,:)=Hfn(n,:)+Hfsegn0;
        end
    end
end

plot3(0,0,0); hold on;
Lmin = min(seglen);
LINKMAX=length(links(:,1));
for i=1:LINKMAX
    n0=links(i,1);
    n1=links(i,2);
    if((n0~=0)&(n1~=0)&(max(max(abs(rn([n0,n1],:))))<=plim))
        lvec = (rn(n1,1:3)-rn(n0,1:3));
        bvec = links(i,3:5);
%         color = 'white';
        color = 'k';
        dline = norm(lvec);
        linedir=lvec./dline;
        costh2=(linedir*bvec')^2/(bvec*bvec');
        sinth2 = 1-costh2;
        if sinth2 < 0.0076%HY20180503: +-5 degrees
%         if sinth2 < 1E-4%HY20180503: +-5 degrees
            color = 'r';%HY20180504: plot pure screw components with color red
        end
        
        PKscaler = 0.;
        HPKscaler = 0.;
        if sinth2 < 0.0076%HY20180503: +-5 degrees
%         if sinth2 < 1E-4%HY20180503: +-5 degrees
            color = 'r';%HY20180504: plot pure screw components with color red
            PKscaler = 0.5E4;
            HPKscaler = 5E4;
        end  
%         if (abs(fn(n1,1))>abs(fn(n1,2))*0.1)|(abs(Hfn(n1,1))>abs(Hfn(n1,2))*0.1)
%             PKscaler = 0.;
%             HPKscaler = 0.;
%             color = 'k';
%         end
%         if norm(Hfn)>norm(fn)*0.1
%             PKscaler = 0.;
%             HPKscaler = 0.;
%             color = 'k';
%         end
%         if rn(n0,1)<300|rn(n1,1)<300|rn(n0,1)>1300|rn(n1,1)>1300|rn(n1,1)<0
% %         if rn(n1,1)<0
%             PKscaler = 0.;
%             HPKscaler = 0.;
%             color = 'k';
%             flag(n0,1) = 1;
%             flag(n1,1) = 1;
%         end
        
%         if rn(n0,1)<300|rn(n0,1)>1300
%             PKscaler = 0.;
%             HPKscaler = 0.;
%             color = 'k';
%             flag(n0,1) = 1;
%         end

        
%         if rn(n0,1)<200|rn(n1,1)<200|rn(n0,1)>2700|rn(n1,1)>2700|rn(n1,1)<0
% %         if rn(n1,1)<0
%             PKscaler = 0.;
%             HPKscaler = 0.;
%             color = 'k';
%             flag(n0,1) = 1;
%             flag(n1,1) = 1;
%         end
%         
%         if rn(n0,1)<300|rn(n0,1)>2700
% %         if rn(n1,1)<0
%             PKscaler = 0.;
%             HPKscaler = 0.;
%             color = 'k';
%             flag(n0,1) = 1;
%         end

        
        
        %filter out "infinity" lines
%         plot3(rn(n0,1),rn(n0,2),rn(n0,3),'ko','LineWidth',3);
%         plot3(rn([n0,n1],1),rn([n0,n1],2),rn([n0,n1],3),'k-','LineWidth',3);
        plot3(rn([n0,n1],1),rn([n0,n1],2),rn([n0,n1],3),color,'LineWidth',4);
%         r0 = rn(n0,1:3);
%         lvec = rn(n1,1:3)-rn(n0,1:3);
%         quiver3(r0(1),r0(2),r0(3),lvec(1),lvec(2),lvec(3),'k-','LineWidth',1); 
        Nint = round(norm(lvec)/(Lmin/3));%HY20180603: interpolation points
%         Nint = 1;%HY20180603: interpolation points
        deltafn = (fn(n1,1:3)-fn(n0,1:3))/Nint;
%         deltaHfn = (Hfn(n1,1:3)-Hfn(n0,1:3))/Nint;
        deltarn = (rn(n1,1:3)-rn(n0,1:3))/Nint;
        for k=0:Nint
            quiver3(rn(n0,1)+k*deltarn(1),rn(n0,2)+k*deltarn(2),rn(n0,3)+k*deltarn(3),...
                fn(n0,1)*PKscaler+k*deltafn(1)*PKscaler,fn(n0,2)*PKscaler+k*deltafn(2)*PKscaler,fn(n0,3)*PKscaler+k*deltafn(3)*PKscaler,0, 'g-','LineWidth',2)
%             quiver3(rn(n0,1)+k*deltarn(1),rn(n0,2)+k*deltarn(2),rn(n0,3)+k*deltarn(3),...
%                 Hfn(n0,1)*HPKscaler+k*deltaHfn(1)*HPKscaler,Hfn(n0,2)*HPKscaler+k*deltaHfn(2)*HPKscaler,Hfn(n0,3)*HPKscaler+k*deltaHfn(3)*HPKscaler,0, 'b','LineWidth',4)
        end
    end 
end
amag = 2.856e-4; 
LENGTH_b = sqrt(3)/2;
mumag = 83E3; % MPa only used for plotting 
fnT = fn';
% HfnT = Hfn';
fnnorm = vecnorm(fnT(1:3,find(flag==0))).*amag^2.*mumag;
% Hfnnorm = vecnorm(HfnT(1:3,find(flag==0)))*amag^2*mumag;
max(fnnorm)
% max(Hfnnorm)


% r = 150; % pixels per inch
% % set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1 1]);
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 1600 900]/r);

plot3([-0, 0]',[-plim, plim]',[0, 0]','r--','LineWidth',3);

hold off
axis equal
grid
% set(gca,'fontsize',20)
xlabel('x / b','FontSize',15);
ylabel('y / b','FontSize',15);
xlim([-plim plim]); ylim([-plim plim]); zlim([-plim plim]);