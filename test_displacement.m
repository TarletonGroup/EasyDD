%test displacement script


%% prismatic loop
% rn = [ 2000 -10^4 10^4 0;
% 2100 -10^4 10^4 0;
% 2100 10^4 10^4 0;
% 2000 10^4 10^4 0];
%clear all
rn = [ 2000 2000 -10^5 0;
2100 2000 -10^5 0;
2100 2100 -10^5 0;
2000 2100 -10^5 0;
2000 2000 10^5 0;
2100 2000 10^5 0;
2100 2100 10^5 0;
2000 2100 10^5 0;
];

% rn = [ 2000 2000 -10^5 0;
%     2000 1999 -10^5 0;
% 2100 2000 -10^5 0;
% 2100 1999 -10^5 0;
% 2100 2100 -10^5 0;
% 2100 2099 -10^5 0;
% 2000 2100 -10^5 0;
% 2000 2099 -10^5 0;
% 2000 2000 10^5 0;
% 2000 1999 10^5 0;
% 2100 2000 10^5 0;
% 2100 1999 10^5 0;
% 2100 2100 10^5 0;
% 2100 2099 10^5 0;
% 2000 2100 10^5 0;
% 2000 2099 10^5 0;
% ];


links=[1 2 0 0 1 0 0 0 ;
2 3 0 0 1 0 0 0;
3 4 0 0 1 0 0 0;
4 1 0 0 1 0 0 0;

5 8 0 0 1 0 0 0;
8 7 0 0 1 0 0 0;
7 6 0 0 1 0 0 0;
6 5 0 0 1 0 0 0;
];

% links = [ 2 3 0 0 1 0 0 0 ;
%     3 11 0 0 1 0 0 0 ;
%     11 10 0 0 1 0 0 0 ;
%     10 2 0 0 1 0 0 0 ;
%     
%     1 9 0 0 1 0 0 0 ;
%     9 16 0 0 1 0 0 0 ;
%     16 8 0 0 1 0 0 0 ;
%     8 1 0 0 1 0 0 0 ;
%     
%     7 15 0 0 1 0 0 0 ;
%     15 14 0 0 1 0 0 0 ;
%     14 6 0 0 1 0 0 0 ;
%     6 7 0 0 1 0 0 0 ;
%     
%     4 5 0 0 1 0 0 0 ;
%     5 13 0 0 1 0 0 0 ;
%     13 12 0 0 1 0 0 0 ;
%     12 4 0 0 1 0 0 0];
 
maxconnections=4;
segments=constructsegmentlist(rn,links);
loop_list = loop_search(rn,links,maxconnections);
MU=1;NU=0.305;

[X,Y] = meshgrid(1902:5:2202,1902:5:2202);
Z=zeros(size(X,1),size(Y,1));

% [ut] = displacement_fivel([0,0,0],segments,NU)

for i = 1:size(X,1)
    for j=1:size(Y,1)
        %for shear loop, loop at where x=0, and y,z vary
        %[ut] = displacement([X(i,j),Y(i,j),Z(i,j)],segments,loop_list,NU);
        [ut] = displacement_fivel([X(i,j),Y(i,j),Z(i,j)],segments,NU);
        utilda(i,j,1:3)=ut;
        %omega(i,j,1:4)=O;
    end
end

surf(X,Y,utilda(:,:,3)); %z-displacement
% hold on
% plot3(rn(1,1),rn(1,2),rn(1,3),'+');
% plot3(rn(2,1),rn(2,2),rn(2,3),'+');
% plot3(rn(3,1),rn(3,2),rn(3,3),'+');
% plot3(rn(4,1),rn(4,2),rn(4,3),'+');
% xlabel('x-direction','FontSize',14);
% ylabel('y-direction','FontSize',14);
% zlabel('z-direction','FontSize',14);
% xlim([2000 2100]);
% ylim([2000 2100]);
% title('prismatic loop (b=001), z-direction displacement','FontSize',14);

%% shear loop
%clear all
rn = [ 10 0 -100 0 ;
10 0 100 0;
10 20 100 0;
10 20 -100 0];

links=[1 2 0 0 1 0 1 0 ;
2 3 0 0 1 0 1 0;
3 4 0 0 1 0 1 0;
4 1 0 0 1 0 1 0];

maxconnections=4;
segments=constructsegmentlist(rn,links);
loop_list = loop_search(rn,links,maxconnections);
MU=1;NU=0.305;

% nodes = [10.4082   10.4082    9.1835
%      9.7011   11.1154    9.1835
%      9.2929   10.7071   10.0000
%      8.8846   10.2989   10.8165
% %      9.5918    9.5918   10.8165
%      10.2989    8.8846   10.8165
%      10.7071    9.2929   10.0000
%      11.1154    9.7011    9.1835];
% 
%  b = [0 0 1];

[X,Y] = meshgrid(0.1:0.25:20.1,0.1:0.25:20.1);
Z=zeros(size(X,1),size(Y,1));

x0 = horzcat( reshape(X,size(X,1)*size(Y,1),1) , reshape(Y,size(X,1)*size(Y,1),1) , reshape(Z,size(X,1)*size(Y,1),1) );

[Ux, Uy, Uz] = displacement_fivel(x0,segments,NU);

ux = reshape(Ux,size(X,1),size(Y,1));
uy = reshape(Uy,size(X,1),size(Y,1));
uz = reshape(Uz,size(X,1),size(Y,1));

uz_theory = (1/(2*pi))*atan(y/x);
surf(X,Y,uz); %z-displacement
hold on
% plot3(nodes(1,1),nodes(1,2),nodes(1,3),'+');
% plot3(nodes(2,1),nodes(2,2),nodes(2,3),'+');
% plot3(nodes(3,1),nodes(3,2),nodes(3,3),'+');
% plot3(nodes(4,1),nodes(4,2),nodes(4,3),'+');
% xlabel('x-direction','FontSize',14);
% ylabel('y-direction','FontSize',14);
% zlabel('z-direction','FontSize',14);
% title('shear loop (b=001), x-direction displacement','FontSize',14);
% plot3(nodes(5,1),nodes(5,2),nodes(5,3),'+');
% plot3(nodes(6,1),nodes(6,2),nodes(6,3),'+');