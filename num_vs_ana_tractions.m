%Check that left face of cantilever has same ftilda_n for FEM and analytical.
%clear all;
% %load a test condition
% setenv('MW_NVCC_PATH','C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v7.5\bin')
%
% mexcuda COMPFLAGS="$COMPFLAGS -arch=sm_50 -use_fast_math" ./nodalForce/nodal_surface_force_linear_rectangle_mex.cu -Dsc
clear all
close all
run ./Inputs/i_prismatic_bcc.m
min_mx = 20;
max_mx = 20;
stp_mx = 10;

fig0 = plotnodes(rn,links,0,vertices)
for mx = min_mx: stp_mx: max_mx
    %% FEM
    segments = constructsegmentlist(rn,links);
    [B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
        Sfront,Sback,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
        w,h,d,my,mz,mel] = finiteElement3D(dx,dy,dz,mx,MU,NU,loading);

    sizeRectangleList = mx*mz;
    sizeSegmentList = size(segments,1);
    
    %% Analytical.
    x1_array = reshape(segments(:,6:8)',sizeSegmentList*3,1);
    x2_array = reshape(segments(:,9:11)',sizeSegmentList*3,1);
    b_array  = reshape(segments(:,3:5)',sizeSegmentList*3,1);  
    
    x3_array = zeros(sizeRectangleList,3);
    x4_array = zeros(sizeRectangleList,3);
    x5_array = zeros(sizeRectangleList,3);
    x6_array = zeros(sizeRectangleList,3);
    
    for i=1:sizeRectangleList
        x3x6 = xnodes(nc(i,[5,6,8,7]),1:3);
        x3_array(i,1:3) = x3x6(1,1:3);
        x4_array(i,1:3) = x3x6(2,1:3);
        x5_array(i,1:3) = x3x6(3,1:3);
        x6_array(i,1:3) = x3x6(4,1:3);
    end
    x3_array = reshape(x3_array',sizeRectangleList*3,1); %x,y,z for each node
    x4_array = reshape(x4_array',sizeRectangleList*3,1);
    x5_array = reshape(x5_array',sizeRectangleList*3,1);
    x6_array = reshape(x6_array',sizeRectangleList*3,1);

     %% Serial
     tic;
      [fx3_array,fx4_array,fx5_array,fx6_array,fxtot_array] = nodal_surface_force_linear_rectangle_mex(x1_array,x2_array,...
        x3_array,x4_array,x5_array,x6_array,...
        b_array,MU,NU,a,sizeRectangleList,sizeSegmentList);
    time_a_lin = toc;
    ftilda_a_lin = reshape(fxtot_array,3,mx*mz)';
    
%     %% Parallel
%     tic;
      %[fx3_array,fx4_array,fx5_array,fx6_array,fxtot_array] = nodal_surface_force_linear_rectangle_mex_cuda(x1_array,x2_array,...
       % x3_array,x4_array,x5_array,x6_array,...
        %b_array,MU,NU,a,sizeRectangleList,sizeSegmentList, 512, 1);
%     time_a_par = toc;
%     ftilda_a_par = reshape(fxtot_array,3,mx*mz)';
    
    %% Numerical
    midpoint_element = zeros(mx*mz,3);
    for p=1:mx*mz
       x3x6 = xnodes(nc(p,[5,6,8,7]),1:3);
       midpoint_element(p,1:3) = mean(x3x6,1);
    end %for
    tic;
    % min(y) xz plane
    gamma = gammat(gammat(:,4) == -1,:); % face normal = [0 -1 0];
    area = zeros(mx*mz,1) + gamma(1,2);
    normal = gamma(1,3:5);
    lseg = size(segments,1);
    lgrid = mx*mz;
    p1x = segments(:,6);
    p1y = segments(:,7);
    p1z = segments(:,8);
    p2x = segments(:,9);
    p2y = segments(:,10);
    p2z = segments(:,11);
    bx = segments(:,3);
    by = segments(:,4);
    bz = segments(:,5);
    x = midpoint_element(:,1);
    y = midpoint_element(:,2);
    z = midpoint_element(:,3);
    tic;
    [sxx, syy, szz, sxy, syz, sxz] = StressDueToSegs(lgrid, lseg,...
                                                  x, y, z,...
                                                  p1x,p1y,p1z,...
                                                  p2x,p2y,p2z,...
                                                  bx,by,bz,...
                                                  a,MU,NU); 
    Tx = sxx*normal(1) + sxy*normal(2) + sxz*normal(3);
    Ty = sxy*normal(1) + syy*normal(2) + syz*normal(3);
    Tz = sxz*normal(1) + syz*normal(2) + szz*normal(3);
    ATx = area.*Tx;
    ATy = area.*Ty;
    ATz = area.*Tz;
    ftilda_n = [ATx,ATy,ATz];
    time_n = toc;
    
    %% Plotting.
    X = reshape(x,mx,mz);
    Y = reshape(y,mx,mz);
    Z = reshape(z,mx,mz);
    
    fig1 = figure;
    rel_err = (ftilda_n - ftilda_a_lin)./ftilda_a_lin;
    rel_err(isnan(rel_err)) = 0;
    
    subplot(3,1,1)
    min_rel_err = min(rel_err(:,1))
    max_rel_err = max(rel_err(:,1))
    mean_rel_err = mean(rel_err(:,1))
    contourf(X,Z,reshape(rel_err(:,1),mx,mz));
    colorvec = [min_rel_err max_rel_err];
    caxis(colorvec);
    colorbar
    title('Relative Error, $\mathbf{F_{x}}$','Interpreter','latex');
    ylabel('$z$','Interpreter','latex');
    xlabel('$x$', 'Interpreter', 'latex');
    
    subplot(3,1,2)
    min_rel_err = min(rel_err(:,2))
    max_rel_err = max(rel_err(:,2))
    mean_rel_err = mean(rel_err(:,2))
    contourf(X,Z,reshape(rel_err(:,2),mx,mz));
    colorvec = [min_rel_err max_rel_err];
    caxis(colorvec);
    colorbar
    title('Relative Error, $\mathbf{F_{y}}$','Interpreter','latex');
    ylabel('$z$','Interpreter','latex');
    xlabel('$x$', 'Interpreter', 'latex');
    
    subplot(3,1,3)
    min_rel_err = min(rel_err(:,3))
    max_rel_err = max(rel_err(:,3))
    mean_rel_err = mean(rel_err(:,3))
    contourf(X,Z,reshape(rel_err(:,3),mx,mz));
    colorvec = [min_rel_err max_rel_err];
    caxis(colorvec);
    colorbar
    title('Relative Error, $\mathbf{F_{z}}$','Interpreter','latex');
    ylabel('$z$','Interpreter','latex');
    xlabel('$x$', 'Interpreter', 'latex');
    
    
    
    h=figure;
    subplot(3,2,1);
    
    contourf(X,Z,reshape(ftilda_n(:,1),mx,mz));
%     colorvec = [-0.06 0.06];
%     caxis(colorvec);
    colorbar
    title(['mx = ',num2str(mx)])
    title('$\mathbf{f}^{\mathrm{N}}_x$', 'Interpreter', 'latex');
    xlabel('$x\, (\mu\mathrm{m})$', 'Interpreter', 'latex');
    subplot(3,2,2);
    contourf(X,Z,reshape(ftilda_a_lin(:,1),mx,mz));
%     caxis(colorvec);
    colorbar
    title('$\mathbf{f}^{\mathrm{A}}_x$', 'Interpreter', 'latex');
    ylabel('$z\, (\mu\mathrm{m})$', 'Interpreter', 'latex')
    xlabel('$x\, (\mu\mathrm{m})$', 'Interpreter', 'latex')
    subplot(3,2,3);
    contourf(X,Z,reshape(ftilda_n(:,2),mx,mz));
%     caxis(colorvec);
    colorbar
    title('$\mathbf{f}^{\mathrm{N}}_y$', 'Interpreter', 'latex');
    xlabel('$x\, (\mu\mathrm{m})$', 'Interpreter', 'latex')
    subplot(3,2,4);
    contourf(X,Z,reshape(ftilda_a_lin(:,2),mx,mz));
%     caxis(colorvec);
    colorbar
    title('$\mathbf{f}^{\mathrm{A}}_y$', 'Interpreter', 'latex');
    ylabel('$z\, (\mu\mathrm{m})$', 'Interpreter', 'latex')
    xlabel('$x\, (\mu\mathrm{m})$', 'Interpreter', 'latex')
    subplot(3,2,5);
    contourf(X,Z,reshape(ftilda_n(:,3),mx,mz));
%     caxis(colorvec);
    colorbar
    title('$\mathbf{f}^{\mathrm{N}}_z$', 'Interpreter', 'latex');
    xlabel('$x\, (\mu\mathrm{m})$', 'Interpreter', 'latex')
    subplot(3,2,6);
    contourf(X,Z,reshape(ftilda_a_lin(:,3),mx,mz));
%     caxis(colorvec);
    colorbar
    title('$\mathbf{f}^{\mathrm{A}}_z$', 'Interpreter', 'latex');
    ylabel('$z\, (\mu\mathrm{m})$', 'Interpreter', 'latex')
    xlabel('$x\, (\mu\mathrm{m})$', 'Interpreter', 'latex')
    saveas(h, sprintf('contour%d',mx), 'epsc');
end %for

for mx = min_mx: stp_mx: max_mx
    %% Use the xy plane for min z as the other benchmark
clear x3_array x4_array x5_array x6_array;
clear fx3_array fx4_array fx5_array fx6_array fxtot_array;
clear ftilda_a midpoint_element rel_err;

    %% use mex version of above
    
    segments = constructsegmentlist(rn,links);
    sizeSegmentList = size(segments,1);
    
    x1_array = reshape(segments(:,6:8)',sizeSegmentList*3,1);
    x2_array = reshape(segments(:,9:11)',sizeSegmentList*3,1);
    b_array  = reshape(segments(:,3:5)',sizeSegmentList*3,1); 

    segments = constructsegmentlist(rn,links);
    %construct finite element arrays
    %construct stiffeness matrix K and pre-compute L,U decompositions.
    disp('Constructing stiffness matrix K and precomputing L,U decompositions. Please wait.'); 
    [B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
        Sfront,Sback,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
        w,h,d,my,mz,mel] = finiteElement3D(dx,dy,dz,mx,MU,NU,loading);
    
    midpoint_element = zeros(mx*my,3);
    sizeRectangleList = mx*my;

x3_array=zeros(sizeRectangleList,3);
x4_array=zeros(sizeRectangleList,3);
x5_array=zeros(sizeRectangleList,3);
x6_array=zeros(sizeRectangleList,3);

idx = 1;
for i = 1: my
    cntr = (i-1)*mx*mz;
    for j = 1: mx
        x3x6 = xnodes(nc(cntr + j,[8, 7, 4, 3]),1:3);
        x3_array(idx,1:3) = x3x6(1,1:3);
        x4_array(idx,1:3) = x3x6(2,1:3);
        x5_array(idx,1:3) = x3x6(3,1:3);
        x6_array(idx,1:3) = x3x6(4,1:3);
        midpoint_element(idx,1:3) = mean(x3x6,1);
        idx = idx + 1;
    end %for 
end %for

x3_array = reshape(x3_array',sizeRectangleList*3,1); %x,y,z for each node
x4_array = reshape(x4_array',sizeRectangleList*3,1);
x5_array = reshape(x5_array',sizeRectangleList*3,1);
x6_array = reshape(x6_array',sizeRectangleList*3,1);

     tic;
      [fx3_array,fx4_array,fx5_array,fx6_array,fxtot_array] = nodal_surface_force_linear_rectangle_mex(x1_array,x2_array,...
        x3_array,x4_array,x5_array,x6_array,...
        b_array,MU,NU,a,sizeRectangleList,sizeSegmentList);
    time_lin = toc;
    ftilda_a_lin = reshape(fxtot_array,3,mx*my)';
    
%     tic;
      %[fx3_array2,fx4_array2,fx5_array2,fx6_array2,fxtot_array2] = nodal_surface_force_linear_rectangle_mex_cuda(x1_array,x2_array,...
       % x3_array2,x4_array2,x5_array2,x6_array2,...
        %b_array,MU,NU,a,sizeRectangleList,sizeSegmentList, 512, 1);
%     time_par = toc;
        
    
    
    
    
    
    %%
    tic;
    %choose front face
    gamma = gammat(gammat(:,5) == 1,:); %face=[0 0 -1];
    area = zeros(mx*my,1) + gamma(1,2);
    normal = gamma(1,3:5);
    lseg = size(segments,1);
    lgrid = mx*my;
    p1x = segments(:,6);
    p1y = segments(:,7);
    p1z = segments(:,8);
    p2x = segments(:,9);
    p2y = segments(:,10);
    p2z = segments(:,11);
    bx = segments(:,3);
    by = segments(:,4);
    bz = segments(:,5);
    % x = xnodes(nodes,1);
    % y = xnodes(nodes,2);
    % z = xnodes(nodes,3);
    x = midpoint_element(:,1);
    y = midpoint_element(:,2);
    z = midpoint_element(:,3);
    [sxx, syy, szz, sxy, syz, sxz] = StressDueToSegs(lgrid, lseg,...
                                                  x, y, z,...
                                                  p1x,p1y,p1z,...
                                                  p2x,p2y,p2z,...
                                                  bx,by,bz,...
                                                  a,MU,NU); 
    Tx = sxx*normal(1) + sxy*normal(2) + sxz*normal(3);
    Ty = sxy*normal(1) + syy*normal(2) + syz*normal(3);
    Tz = sxz*normal(1) + syz*normal(2) + szz*normal(3);
    ATx = area.*Tx;
    ATy = area.*Ty;
    ATz = area.*Tz;
    ftilda_n = [ATx,ATy,ATz];

    X=reshape(x,mx,my);
    Y=reshape(y,mx,my);
    Z=reshape(z,mx,my);
    toc;
    
    fig1 = figure;
    rel_err = (ftilda_n - ftilda_a_lin)./ftilda_a_lin;
    rel_err(isnan(rel_err)) = 0;
    
    subplot(3,1,1)
    min_rel_err = min(rel_err(:,1))
    max_rel_err = max(rel_err(:,1))
    mean_rel_err = mean(rel_err(:,1))
    contourf(X,Y,reshape(rel_err(:,1),mx,my));
    colorvec = [min_rel_err max_rel_err];
    caxis(colorvec);
    colorbar
    title('Relative Error, $\mathbf{F_{x}}$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    xlabel('$x$', 'Interpreter', 'latex');
    
    subplot(3,1,2)
    min_rel_err = min(rel_err(:,2))
    max_rel_err = max(rel_err(:,2))
    mean_rel_err = mean(rel_err(:,2))
    contourf(X,Y,reshape(rel_err(:,2),mx,my));
    colorvec = [min_rel_err max_rel_err];
    caxis(colorvec);
    colorbar
    title('Relative Error, $\mathbf{F_{y}}$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    xlabel('$x$', 'Interpreter', 'latex');
    
    subplot(3,1,3)
    min_rel_err = min(rel_err(:,3))
    max_rel_err = max(rel_err(:,3))
    mean_rel_err = mean(rel_err(:,3))
    contourf(X,Y,reshape(rel_err(:,3),mx,my));
    colorvec = [min_rel_err max_rel_err];
    caxis(colorvec);
    colorbar
    title('Relative Error, $\mathbf{F_{z}}$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    xlabel('$x$', 'Interpreter', 'latex');
    
    h=figure;
    subplot(3,2,1);
    
    contourf(X,Y,reshape(ftilda_n(:,1),mx,my));
%     colorvec = [-0.06 0.06];
%     caxis(colorvec);
    colorbar
    title(['mx = ',num2str(mx)])
    title('$\mathbf{f}^{\mathrm{N}}_x$', 'Interpreter', 'latex');
    xlabel('$x\, (\mu\mathrm{m})$', 'Interpreter', 'latex');
    subplot(3,2,2);
    contourf(X,Y,reshape(ftilda_a_lin(:,1),mx,my));
%     caxis(colorvec);
    colorbar
    title('$\mathbf{f}^{\mathrm{A}}_x$', 'Interpreter', 'latex');
    ylabel('$z\, (\mu\mathrm{m})$', 'Interpreter', 'latex')
    xlabel('$x\, (\mu\mathrm{m})$', 'Interpreter', 'latex')
    subplot(3,2,3);
    contourf(X,Y,reshape(ftilda_n(:,2),mx,my));
%     caxis(colorvec);
    colorbar
    title('$\mathbf{f}^{\mathrm{N}}_y$', 'Interpreter', 'latex');
    xlabel('$x\, (\mu\mathrm{m})$', 'Interpreter', 'latex')
    subplot(3,2,4);
    contourf(X,Y,reshape(ftilda_a_lin(:,2),mx,my));
%     caxis(colorvec);
    colorbar
    title('$\mathbf{f}^{\mathrm{A}}_y$', 'Interpreter', 'latex');
    ylabel('$z\, (\mu\mathrm{m})$', 'Interpreter', 'latex')
    xlabel('$x\, (\mu\mathrm{m})$', 'Interpreter', 'latex')
    subplot(3,2,5);
    contourf(X,Y,reshape(ftilda_n(:,3),mx,my));
%     caxis(colorvec);
    colorbar
    title('$\mathbf{f}^{\mathrm{N}}_z$', 'Interpreter', 'latex');
    xlabel('$x\, (\mu\mathrm{m})$', 'Interpreter', 'latex')
    subplot(3,2,6);
    contourf(X,Y,reshape(ftilda_a_lin(:,3),mx,my));
%     caxis(colorvec);
    colorbar
    title('$\mathbf{f}^{\mathrm{A}}_z$', 'Interpreter', 'latex');
    ylabel('$z\, (\mu\mathrm{m})$', 'Interpreter', 'latex')
    xlabel('$x\, (\mu\mathrm{m})$', 'Interpreter', 'latex')
    saveas(h, sprintf('contour%d',mx), 'epsc');
end %for

