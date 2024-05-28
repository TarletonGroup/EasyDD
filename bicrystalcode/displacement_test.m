

segments = constructsegmentlist(rn,links);

x0 = [linspace(0,dx)', 2.5e4*ones(100,1), 4.5e4*ones(100,1)];

tic;
[Ux_f,Uy_f,Uz_f]  = displacement_fivel(x0,segments,NU);
toc;

point_array_length = size(x0,1);
segments_array_length = size(segments,1);
tic;
[Ux,Uy,Uz] = UtildaMex(x0(:,1),x0(:,2),x0(:,3),... %coordinates
                       segments(:,3), segments(:,4), segments(:,5),... %burgers vector
                       segments(:,6), segments(:,7), segments(:,8),... %start node segs
                       segments(:,9), segments(:,10), segments(:,11),... %end node segs
                       segments(:,12), segments(:,13), segments(:,14),... %slip plane
                       NU,point_array_length,segments_array_length);
utilda2 = horzcat(Ux,Uy,Uz);
toc;

norm([Ux-Ux_f,Uy-Uy_f,Uz-Uz_f]);

clear utilda3
for i =1:size(x0,1)
    
uspf = displacement_spf(x0(i,:),rn(:,1:3),segments(1,3:5),NU);
utilda3(i,1:3) = uspf;
end

