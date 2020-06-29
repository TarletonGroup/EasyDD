%Check that left face of cantilever has same ftilda_n for FEM and analytical.
%clear all;
% %load a test condition
% setenv('MW_NVCC_PATH','C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.0\bin')
% mexcuda COMPFLAGS="$COMPFLAGS -arch=sm_70 -use_fast_math" ./nodalForce/nodal_surface_force_linear_rectangle_mex_cuda.cu
% The above should work, but doesn't this is a known matlab issue. I
% fixed it by following the instructions found here:
% https://devtalk.nvidia.com/default/topic/1027876/cuda-programming-and-performance/why-does-atomicadd-not-work-with-doubles-as-input-/post/5228325/#5228325
% clear all
% close all
run ./Inputs/i_prismatic_bcc.m
min_mx = 30;
max_mx = 30;
stp_mx = 10;

% fig0 = plotnodes(rn,links,0,vertices)
for mx = min_mx: stp_mx: max_mx
    %close all
    clear x3_array x4_array x5_array x6_array;
    clear fx3_array fx4_array fx5_array fx6_array fxtot_array;
    clear ftilda_a midpoint_element rel_err;
    %% FEM
    segments = constructsegmentlist(rn,links);
    [B,xnodes,mno,nc,n,D,kg,K,L,U,Sleft,Sright,Stop,Sbot,...
        Sfront,Sback,gammat,gammau,gammaMixed,fixedDofs,freeDofs,...
        w,h,d,my,mz,mel] = finiteElement3D(dx,dy,dz,mx,MU,NU,loading);

    sizeSegmentList = size(segments,1);
    x1_array = reshape(segments(:,6:8)',sizeSegmentList*3,1);
    x2_array = reshape(segments(:,9:11)',sizeSegmentList*3,1);
    b_array  = reshape(segments(:,3:5)',sizeSegmentList*3,1);
    %%
    for face = 1: 1: 3
        % Analytical.
        [dim, x3_array, x4_array, ...
         x5_array, x6_array, midpoint_element] = ...
                                get_face_nodes(face, mx, my, mz, xnodes, nc);
          sizeRectangleList = dim(1);

        %% Serial
        tic;
        [fx3_array, fx4_array, ...
         fx5_array, fx6_array, ...
         fxtot_array] = ...
           nodal_surface_force_linear_rectangle_mex(x1_array, x2_array, ...
              x3_array, x4_array, x5_array, x6_array, b_array, MU, NU, a, ...
              sizeRectangleList, sizeSegmentList, 0);
        time_a_lin = toc;
        ftilda_a_lin = reshape(fxtot_array, 3, dim(1))';

        % Parallel
        tic;
          [fx3_array,fx4_array,fx5_array,fx6_array,fxtot_array] = nodal_surface_force_linear_rectangle_mex_cuda(x1_array,x2_array,...
           x3_array,x4_array,x5_array,x6_array,...
            b_array,MU,NU,a,sizeRectangleList,sizeSegmentList, 128, 1, 0);
        time_a_par = toc;
        ftilda_a_par = reshape(fxtot_array,3,dim(1))';
%         ftilda_mat = NodalSurfForceLinearRectangle2(                    ...
%                                 x1_array, x2_array,...
%                                 x3_array, x4_array,...
%                                 x5_array, x6_array,...
%                                 b_array, MU,NU,a);

        maxval = max(ftilda_a_lin./ftilda_a_par)
        minval = min(ftilda_a_lin./ftilda_a_par)
        tdiff = time_a_lin/time_a_par

        % Numerical
        [ftilda_n, x, y, z] = f_dln_num(face, dim, gammat, segments, midpoint_element, a, MU, NU);
%         max(ftilda_n - ftilda_a_par)
%         min(ftilda_n - ftilda_a_par)
%         mean(ftilda_n - ftilda_a_par)

        % Plotting
%         save(sprintf('mx=%d_face=%d', mx, face), 'ftilda_a_lin', 'ftilda_a_par', 'ftilda_n', 'dim', 'face', 'midpoint_element')
%         plot_figs(face, dim, ftilda_n, ftilda_a_lin, x, y, z);
    end %for

end %for

%% Plot
min_face = 1;
stp_face = 1;
max_face = 3;
n_plots = (max_mx - min_mx)/stp_mx;
rms_rel_err_val = zeros(max_face, (max_mx - min_mx + stp_mx)/stp_mx, 3);
cntr = 1;
for face = min_face: stp_face: max_face
    for mx = 100: stp_mx: max_mx
        load(sprintf('mx=%d_face=%d', mx, face));
        rms_rel_err_val(face, cntr, :) = ...
                        plot_figs(face, dim, ftilda_n, ftilda_a_lin, ...
                        midpoint_element(:, 1), midpoint_element(:, 2),...
                        midpoint_element(:, 3));
        cntr = cntr + 1;
    end %for
end %for
x = min_mx: stp_mx: max_mx;

n_plots = n_plots+1;
for face = min_face: stp_face: max_face
    init = 1 + (face-min_face)*n_plots;
    fin = init + n_plots - 1;
    Y = [rms_rel_err_val(face, : , 1); rms_rel_err_val(face, : , 2); rms_rel_err_val(face, : , 3)];
    rms_fig = figure;
    plot(x, Y)
    legend({'$x$', '$y$', '$z$'}, 'Interpreter', 'latex');
end %for
%%
function [x3, x4, x5, x6, xmid] = establish_se(length)
    x3   = zeros(length, 3);
    x4   = zeros(length, 3);
    x5   = zeros(length, 3);
    x6   = zeros(length, 3);
    xmid = zeros(length, 3);
end %function

function [x3, x4, x5, x6] = reshape_se(x3, x4, x5, x6, length)
    x3 = reshape(x3', length*3, 1);
    x4 = reshape(x4', length*3, 1);
    x5 = reshape(x5', length*3, 1);
    x6 = reshape(x6', length*3, 1);
end %function

function [dim, x3, x4, x5, x6, xmid] = get_face_nodes(face, mx, my, mz, xnodes, nc)
    if face == 1
        dim = [mx*mz; mx; mz];
        [x3, x4, x5, x6, xmid] = establish_se(dim(1));
        nodes = [5, 6, 8, 7];
        for i=1: dim(1)
            x3x6 = xnodes(nc(i, nodes), 1:3);
            x3  (i, 1:3) = x3x6(1   , 1:3);
            x4  (i, 1:3) = x3x6(2   , 1:3);
            x5  (i, 1:3) = x3x6(3   , 1:3);
            x6  (i, 1:3) = x3x6(4   , 1:3);
            xmid(i, 1:3) = mean(x3x6,   1);
        end
    elseif face == 2
        dim = [mx*my; mx; my];
        [x3, x4, x5, x6, xmid] = establish_se(dim(1));
        nodes = [8, 7, 4, 3];
        idx = 1;
        for i = 1: my
            cntr = (i-1)*mx*mz;
            for j = 1: mx
                x3x6 = xnodes(nc(cntr + j, nodes),1:3);
                x3  (idx, 1:3) = x3x6(1   , 1:3);
                x4  (idx, 1:3) = x3x6(2   , 1:3);
                x5  (idx, 1:3) = x3x6(3   , 1:3);
                x6  (idx, 1:3) = x3x6(4   , 1:3);
                xmid(idx, 1:3) = mean(x3x6,   1);
                idx = idx + 1;
            end %for
        end %for
    elseif face == 3
        dim = [my*mz; my; mz];
        [x3, x4, x5, x6, xmid] = establish_se(dim(1));
        nodes = [6, 2, 7, 3];
        idx = 1;
        for i = 1: my
            cntr = (i-1)*mx*mz;
            for j = 1: mz
                x3x6 = xnodes(nc(cntr + j*mx, nodes),1:3);
                x3  (idx, 1:3) = x3x6(1   ,1:3);
                x4  (idx, 1:3) = x3x6(2   ,1:3);
                x5  (idx, 1:3) = x3x6(3   ,1:3);
                x6  (idx, 1:3) = x3x6(4   ,1:3);
                xmid(idx, 1:3) = mean(x3x6,  1);
                idx = idx + 1;
            end %for
        end %for
    end %if
    [x3, x4, x5, x6] = reshape_se(x3, x4, x5, x6, dim(1));
end %function

function [f_dln, x, y ,z] = f_dln_num(face, dim, gammat, segments, xmid, a, MU, NU)
    if face == 1
        idx = 4;
        val = -1;
    elseif face == 2
        idx = 5;
        val = 1;
    elseif face == 3
        idx = 3;
        val = 1;
    end %if
    gamma = gammat(gammat(:, idx) == val,:);
    area = zeros(dim(1),1) + gamma(1,2);
    normal = gamma(1,3:5);
    lseg = size(segments,1);
    lgrid = dim(1);
    p1x = segments(:,6);
    p1y = segments(:,7);
    p1z = segments(:,8);
    p2x = segments(:,9);
    p2y = segments(:,10);
    p2z = segments(:,11);
    bx = segments(:,3);
    by = segments(:,4);
    bz = segments(:,5);
    x = xmid(:,1);
    y = xmid(:,2);
    z = xmid(:,3);
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
    f_dln = [ATx,ATy,ATz];
end %function

function rms_rel_err = plot_figs(face, dim, f_dln_n, f_dln_a, x, y, z)

    if face == 1
        axis1 = reshape(x, dim(2), dim(3));
        axis2 = reshape(z, dim(2), dim(3));
        str_label  = ['x'; 'z'];
    elseif face == 2
        axis1 = reshape(x, dim(2), dim(3));
        axis2 = reshape(y, dim(2), dim(3));
        str_label  = ['x'; 'y'];
    elseif face == 3
        axis1 = reshape(y, dim(2), dim(3));
        axis2 = reshape(z, dim(2), dim(3));
        str_label  = ['y'; 'z'];
    end %if

    str_dim = ['x'; 'y'; 'z'];

    rel_err = (f_dln_n - f_dln_a)./f_dln_a;
    rel_err(isnan(rel_err)) = 0;

    error_cntr = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
    for i = 1: 3
        subplot(3, 1, i);

        contourf(axis1, axis2, reshape(rel_err(:, i), dim(2), dim(3)));

        min_rel_err  = min (rel_err(:, i));
        max_rel_err  = max (rel_err(:, i));
        colorvec = [min_rel_err max_rel_err];
        caxis(colorvec);
        colorbar

        title(sprintf('Relative Error, $F_{%s}$', str_dim(i)),'Interpreter','latex');
        ylabel(sprintf('$%s$', str_label(2)), 'Interpreter', 'latex');
        xlabel(sprintf('$%s$', str_label(1)), 'Interpreter', 'latex');
        axis equal
    end %for

    f_dln_n_fig = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
    for i = 1: 3
        subplot(3, 1, i);

        contourf(axis1, axis2, reshape(f_dln_n(:, i), dim(2), dim(3)));
        mean_rel_err = mean(abs(rel_err(:, i)));
        cvec = [min(f_dln_a(:, i))*(mean_rel_err + 1) max(f_dln_a(:, i))*(mean_rel_err + 1)];
        colorvec = [min(cvec) max(cvec)];%[min(f_dln_a(:, i))*(mean_rel_err + 1) max(f_dln_a(:, i))*(mean_rel_err + 1)];
        caxis(colorvec);
        colorbar
        title(sprintf('$F^{A}_{%s}$', str_dim(i)),'Interpreter','latex');
        ylabel(sprintf('$%s$', str_label(2)), 'Interpreter', 'latex');
        xlabel(sprintf('$%s$', str_label(1)), 'Interpreter', 'latex');
        axis equal
    end %for

    f_dln_a_fig = figure('DefaultAxesPosition', [0.1, 0.1, 0.8, 0.8]);
    for i = 1: 3
        subplot(3, 1, i);
        contourf(axis1, axis2, reshape(f_dln_a(:, i), dim(2), dim(3)));
        mean_rel_err = mean(abs(rel_err(:, i)));
        cvec = [min(f_dln_a(:, i))*(mean_rel_err + 1) max(f_dln_a(:, i))*(mean_rel_err + 1)];
        colorvec = [min(cvec) max(cvec)];
        caxis(colorvec);
        colorbar
        title(sprintf('$F^{A}_{%s}$', str_dim(i)),'Interpreter','latex');
        ylabel(sprintf('$%s$', str_label(2)), 'Interpreter', 'latex');
        xlabel(sprintf('$%s$', str_label(1)), 'Interpreter', 'latex');
        axis equal
    end %for
    rms_rel_err = rms(rel_err);
end %function
