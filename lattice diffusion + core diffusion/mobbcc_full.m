% moified by gaoyuan 09-09-18

function [vn,vn_climb,fn] = mobbcc_full(fseg,rn,links,connectivity,xnodes,mno,wx,wy,wz,mx,my,mz,mel,nc,ke,Dv,nodelist1,conlist1)

global flag_c flag_b doclimb doglide

L0 = size(nodelist1,1);
L1 = size(rn,1);


vn_c=zeros(L1,3);
fn_c=zeros(L1,3);

if doclimb
    
    %contribution from core diffusion
    if flag_c
        [vn_core,fn_core] = mobbcc_core(fseg,rn,links,connectivity,nodelist1,conlist1);
    else
        vn_core=zeros(L1,3);
        fn_core=vn_core;
    end
    
    %contribution from bulk diffusion
    if flag_b
        [vn_bulk,fn_bulk] = mobbcc_bulk(fseg,rn,links,connectivity,xnodes,mno,wx,wy,wz,mx,my,mz,mel,nc,ke,Dv,nodelist1,conlist1);
    else
        vn_bulk=zeros(L1,3);
        fn_bulk=vn_bulk;
    end
    
    vn_c = vn_core*flag_c + vn_bulk*flag_b;
    
    if flag_c
        fn_c = fn_core;
    else
        fn_c=fn_bulk;
    end
end


vn_g=zeros(L1,3);
fn_g=zeros(L1,3);

if doglide
    [vn_g,fn_g] = mobbcc_glide(fseg,rn,links,connectivity,nodelist1,conlist1);
end


% for the remesh (L0~=0)
if ~L0
    
    if doglide
        fn=fn_g;
        vn=vn_g;
    else
        fn=fn_c;
        vn_climb=vn_c;
        vn=vn_climb;
    end
    
else
    vn_climb=zeros(L0,3);
    vn=zeros(L0,3);
    fn=zeros(L0,3);
    for i=1:L0
        vn_climb(i,:)=vn_c(nodelist1(i),1:3);
        
        if doglide
            fn(i,:)=fn(i,1:3);
            vn(i,:)=vn_g(i,:);
        else
            vn(i,:)=vn_climb(i,:);
            fn(i,:)=fn_c(nodelist1(i),1:3);
        end
    end
end

end