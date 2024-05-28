

function [vn_c,links,fn_c] = mobbcc_linear_full_climb(fseg,rn,links,connectivity,nodelist1,conlist1)

% vn  total velocity
% vn_c climb velocity
% vn_g glide velocity
% fn1 total nodal force

%mobility law function (model: FCC1)
%FCC1: velocity linear to force but remain orthogonal to glide plane
% if nodelist and conlist are empty then vn is the length of rn otherwise
% vn is the length of nodelist

global do_inclusion 

% L0 = size(nodelist1,1);

L1 = size(rn,1);    


nearinc_flag=zeros(L1,1);  % flag for the climbing nodes around the inclusions, flag=1 for nodes on the inclusion and climb toward the center

[v_c,fn_c]=mobbcc_linear_full_climb_test(fseg,rn,links,connectivity);

if do_inclusion ==1
    [nearinc_flag]=inc_velocity_climb(v_c,rn,links,nearinc_flag);
end

% 
vn_c=v_c.*((nearinc_flag(:,1)==0)*[1 1 1]);

% [v_c,fn_c]=mobbcc_linear_full_climb_test(fseg,rn,links,connectivity,nearinc_flag);
rnmax=size(rn,2);
for i=1:size(rn,1)
    if rn(i,rnmax)~=0
        v_c(i,1:3)=0;
    end
end
end





