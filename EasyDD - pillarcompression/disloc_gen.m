%%%%%
%
% Simulation of 20-50 hexagonal <111> loops and square <001> loops in a
% thin film of b.c.c. tungsten
% Lattice vector of tungsten is 0.314 nm. Use this to convert nm into
% lattice vector units for use in DDLab.
%
%%%%%

function [rn,links] = disloc_gen(n_loops,randSizeMean, randSizeDev,Lx_box,Ly_box,thickness)
%Simulation box

bW = 0.314; %conversion nanometers to a

%Loop characteristics
%n_loops = 1000;
%randSizeMean = 3/bW;
%randSizeDev = 0.5/bW;
%Lx_box = 1000/bW;
%Ly_box = 1000/bW;
%thickness = 25/bW; %+/-25nm (total 50nm)
thickness_buffer=thickness;
%-randSizeMean;
Lx = Lx_box;
%+ 2*thickness;
Ly = Ly_box;
%+ 2*thickness;
%randSizeMean;

%% Generate locations for loops
for i = 1:n_loops
    
%generate random loop centerpoint position
loop_clxy = thickness + (Lx_box).*rand(1,2); 
loop_clz = thickness_buffer + (-thickness_buffer-thickness_buffer).*rand(1,1);
loop_loc(i,:) = horzcat(loop_clxy,loop_clz);   

%check this position will not cause intersection with other loops
if i>1
    
    n_loop_exist = size(loop_loc,1);
    current = repmat(loop_loc(i,:),n_loop_exist-1,1);
    distance = abs(current - loop_loc(1:(i-1),:));
    min_distance = min(sqrt(distance(:,1).^2 + distance(:,2).^2 + distance(:,3).^2));
    
    while min_distance < 2*randSizeMean
        
        %re-generate random loop centerpoint position
        loop_clxy = thickness + (Lx_box).*rand(1,2); 
        loop_clz = thickness_buffer + (-thickness_buffer-thickness_buffer).*rand(1,1);
        loop_loc(i,:) = horzcat(loop_clxy,loop_clz);   
  
        current = repmat(loop_loc(i,:),n_loop_exist-1,1);
        distance = current - loop_loc(1:(i-1),:);
        min_distance = min(sqrt(distance(:,1).^2 + distance(:,2).^2 + distance(:,3).^2));

    end

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From X.Yi (2012) there are 25% <001> and 75% <111> at 500C 1.2dpa W     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_001 = 0;
n_111 = 0;

prob=rand(n_loops);
for p=1:size(prob,1)
    if prob(p) <= 8/14
       n_111 = n_111 + 1;
    else
       n_001 = n_001 + 1;
    end
end

n1 = zeros(n_001,3);
n2 = zeros(n_111,3);

rn = zeros(4*n_001+6*n_111,4);
links = zeros(4*n_001+6*n_111,8);

%assign a burgers vector <100>, each combination equally probable.
for p = 1:n_001
    
    randnum = rand(1);
    third = 1/3;
    
    if randnum <= third;
    
    n1(p,:) = [0,0,1];    
        
    elseif randnum <= 2*third && randnum > third
    
    n1(p,:) = [0,1,0];
                 
    elseif randnum <= 3*third && randnum > 2*third
        
    n1(p,:) = [1,0,0];
       
    end
    
end

%assign a burgers vector <111>, each combination equally probable.
for p = 1:n_111
    
    randnum = rand(1);
    eighth = 1/8;
    
    if randnum <= eighth;
    
    n2(p,:) = [1,1,1];    
        
    elseif randnum <= 2*eighth && randnum > eighth
    
    n2(p,:) = [-1,-1,-1];
        
    elseif randnum <= 3*eighth && randnum > 2*eighth
    
    n2(p,:) = [-1,1,1];
       
    elseif randnum <= 4*eighth && randnum > 3*eighth
    
    n2(p,:) = [1,-1,-1];
               
    elseif randnum <= 5*eighth && randnum > 4*eighth
        
    n2(p,:) = [1,-1,1];
       
    elseif randnum <= 6*eighth && randnum > 5*eighth    
    
    n2(p,:) = [-1,1,-1];

    elseif randnum <= 7*eighth && randnum > 6*eighth    
    
    n2(p,:) = [1,1,-1];

    elseif randnum <= 8*eighth && randnum > 7*eighth    
    
    n2(p,:) = [-1,-1,1];

    end
    
end
b = vertcat(n1,n2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now that we've assigned the centerpoint positions for all the loops and %
% their respective burgers vectors, we need to generate the segment       %
% coordinates for use in DDLab. Loop plane perp. to their b-vector.       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Square loop coordinates. Plane normal is 001

sides_square = [ 0,  1, 1;
                 0,  1, -1;
                 0, -1, -1;
                 0, -1, 1];

%Hex loop coordinates. plane normal is 001
fac=sqrt(3)*0.5;
sides_hex = [1,  0,   0;
             0.5,fac, 0;
            -0.5,fac, 0;
            -1,  0,   0;
            -0.5,-fac,0;
             0.5,-fac,0];

%Rotation and translation  
counter=1;
for c=1:n_001
 
v1 = [0,0,1];
v2 = b(c,:);
 
if isequal(abs(v2),[0,0,1])
    
    R = [0 0 1 ; 0 1 0 ; -1 0 0];
    
elseif isequal(abs(v2),[0,1,0])
    
    R = [0 -1 0 ; 1 0 0 ; 0 0 1];
    
elseif isequal(abs(v2),[1,0,0])
    
    R = [1 0 0 ; 0 1 0 ; 0 0 1];
    
end


for z=1:4
    sides_square_trans(z,:) = R*sides_square(z,:)';
end

    coord_square_loops(:,:,counter) = round((randSizeDev*randn)+randSizeMean)*sides_square_trans(:,:) + repmat(loop_loc(c,:),4,1);

    counter=counter+1;
end

counter=1;
for c=(n_001+1):n_loops

v1 = b(c,:);
v2 = [0,0,1];

line=cross(v1,v2);
angle=-acosd(dot(v1,v2)/(norm(v1)*norm(v2)));

for z=1:6

    [sides_hex_trans, R, t] = AxelRot(sides_hex', angle, line, 0);
    
end

coord_hex_loops(:,:,counter) = round((randSizeDev*randn)+randSizeMean)*sides_hex_trans(:,:)' + repmat(loop_loc(c,:),6,1);

counter=counter+1;

end

%
% So now that we have the locations of all the points of all the loops we
% and the burgers vectors we need to from links for DDLab input.
% This has form [x1 y1 z1 0].
%

linksCounter=1;
for j=1:n_001
    
    rn(linksCounter,1:3) = coord_square_loops(1,:,j);
    rn(linksCounter+1,1:3) = coord_square_loops(2,:,j);
    rn(linksCounter+2,1:3) = coord_square_loops(3,:,j);
    rn(linksCounter+3,1:3) = coord_square_loops(4,:,j);
    
    links(linksCounter,1:5) = [linksCounter, linksCounter+1 , b(j,:)];
    links(linksCounter+1,1:5) = [linksCounter+1, linksCounter+2, b(j,:)];
    links(linksCounter+2,1:5) = [linksCounter+2, linksCounter+3, b(j,:)];
    links(linksCounter+3,1:5) = [linksCounter+3, linksCounter, b(j,:)];
    
    %calculate glide planes
    for m=1:4
    line = rn( links(linksCounter+m-1,1) , 1:3) - rn( links(linksCounter+m-1,2) , 1:3);
    plane_n = cross(b(j,:),line);
    links(linksCounter+m-1,6:8) = plane_n;
    end

    linksCounter = linksCounter + 4;
end

for j=(n_001+1):n_loops
    
    rn(linksCounter,1:3) = coord_hex_loops(1,:,j-n_001);
    rn(linksCounter+1,1:3) = coord_hex_loops(2,:,j-n_001);
    rn(linksCounter+2,1:3) = coord_hex_loops(3,:,j-n_001);
    rn(linksCounter+3,1:3) = coord_hex_loops(4,:,j-n_001);
    rn(linksCounter+4,1:3) = coord_hex_loops(5,:,j-n_001);
    rn(linksCounter+5,1:3) = coord_hex_loops(6,:,j-n_001);
    
    links(linksCounter,1:5) = [linksCounter, linksCounter+1 , b(j,:)];
    links(linksCounter+1,1:5) = [linksCounter+1, linksCounter+2, b(j,:)];
    links(linksCounter+2,1:5) = [linksCounter+2, linksCounter+3, b(j,:)];
    links(linksCounter+3,1:5) = [linksCounter+3, linksCounter+4, b(j,:)];
    links(linksCounter+4,1:5) = [linksCounter+4, linksCounter+5, b(j,:)];
    links(linksCounter+5,1:5) = [linksCounter+5, linksCounter, b(j,:)];
    
    %calculate glide planes and update links
    for m=1:6
    line = rn( links(linksCounter+m-1,1) , 1:3) - rn( links(linksCounter+m-1,2) , 1:3);
    plane_n = cross(b(j,:),line);
    links(linksCounter+m-1,6:8) = plane_n;
    end
    
   
    linksCounter = linksCounter + 6;   
end

end






