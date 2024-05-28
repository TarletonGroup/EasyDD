function[RR0,QQ0]=inicon_precipitate(SimBox,volume_frac)

global amag

l_surf = 200*1e-3/amag;   % generate precipitates in a smaller box than the simBox to avoid intersect with surf
crit_dis =10*1e-3/amag; % minmum distance between two precipitates
R=100*1e-3/amag;

% Simbox=SimBox - l_surf;

V=SimBox(1)*SimBox(2)*SimBox(3);
n_sphere = V*volume_frac/(pi*R^3);
n_sphere = floor(n_sphere); % number of sphere

RR0=ones(1,n_sphere)*R;
QQ0=zeros(n_sphere,3);

num_pre=2;
QQ0(1,1)= 0.5*SimBox(1)+(rand-0.5)*(SimBox(1)-l_surf);
QQ0(1,2)=0.5*SimBox(2)+(rand-0.5)*(SimBox(2)-l_surf);
QQ0(1,3)=0.5*SimBox(3)+(rand-0.5)*(SimBox(3)-l_surf);
flag=0;

while num_pre <= n_sphere
 
    QQ0(num_pre,1)= 0.5*SimBox(1)+(rand-0.5)*(SimBox(1)-l_surf);
    QQ0(num_pre,2)= 0.5*SimBox(2)+(rand-0.5)*(SimBox(2)-l_surf);
    QQ0(num_pre,3)= 0.5*SimBox(3)+(rand-0.5)*(SimBox(3)-l_surf);
    for i=1:num_pre-1
        dis=QQ0(num_pre,1:3)-QQ0(i,:);
        dis=norm(dis);
        if (dis-2*R) < crit_dis
            flag=0;
            break
        else
            flag=1;
        end
    end
    if flag==1
        num_pre=num_pre+1;
    end
                
end


fid =   fopen('.\input\initial_structure_precipitate','w');
l1=size(RR0,2);
fwrite(fid,l1,'integer*4');
fwrite(fid,RR0,'real*8');
l2=size(QQ0,1);
fwrite(fid,l2,'integer*4');
fwrite(fid,QQ0,'real*8');
fclose(fid);
end