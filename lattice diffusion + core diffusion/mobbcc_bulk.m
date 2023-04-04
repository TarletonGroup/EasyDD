% modified by gaoyuan 09-09-18

function [vn_b,fn_c] = mobbcc_bulk(fseg,rn,links,connectivity,xnodes,mno,wx,wy,wz,mx,my,mz,mel,nc,ke,Dv,nodelist1,conlist1)


global bur


L1 = size(rn,1); % when mobility law is call from remesh fucnction, L1~=L0;
nodelist = linspace(1,L1,L1)';
[L2,L3] = size(connectivity);
conlist = zeros(L2,(L3-1)/2+1);
conlist(:,1) = connectivity(:,1);
for i=1:L2
    connumb = conlist(i,1);
    conlist(i,2:connumb+1) = linspace(1,connumb,connumb);
end

%% dislocation climb dominated by bulk diffusion
nseg=size(links,1);

v_b=zeros(nseg,3);  % segmental climb velocity
vb=zeros(nseg,1);  % magnitude of v_b

fn=zeros(L1,3); % nodal total force
vn_b=zeros(L1,3); %nodal climb velocity
l=zeros(nseg,1);  %seg length

%% calculate fnc: nodal climb force
fnc=zeros(L1,1);
fn_c=zeros(L1,3);
normc=zeros(L1,3);

for n=1:L1
    n0=nodelist(n);
    numNbrs=conlist(n,1);
    rt=zeros(numNbrs,3);
    
    for i=1:numNbrs
        ii=conlist(n,i+1);
        linkid=connectivity(n0,2*ii);
        posinlink=connectivity(n0,2*ii+1);
        n1=links(linkid,3-posinlink);
        rt=rn(n1,1:3)-rn(n0,1:3);
        
        l(linkid)=norm(rt);
        
        normplane= links(linkid,6:8);               
        % nc(n,:)=nc(n,:)+ normplane*l(linkid);
        normc(n,:)=normc(n,:)+ normplane;
        %             L =L+ l(linkid);        
        if l(linkid)>0
            fsegn0 =fseg(linkid,3*(posinlink-1)+[1:3])/l(linkid);
%             fsegn0 =fseg(linkid,3*(posinlink-1)+[1:3]);
            fn(n,:) = fn(n,:) + fsegn0;
        end
        
    end
    normc(n,:)=normc(n,:)/norm(normc(n,:));
    fnc(n,:)=dot(fn(n,:),normc(n,:));
    fn_c=fnc.*normc;
end

%% calculate k maatrix around dislocation core
[miu,miu_c,Kvc,KEC,KC,gammaMiu_el] = FEMcoupler_new(rn,links,fn_c,normc,xnodes,mno,wx,wy,wz,mx,my,mz,mel,nc,ke,Dv,nodelist,conlist,connectivity);


%---calculate {vc}
%[KEC]'*{miu}+[KC]*{miu_c}=[Kvc]*{vc}
P=KEC'*miu+KC*miu_c;
vc=Kvc\P;
vc=vc/(bur*1e-9)^2;

vn_b=vc.*normc;