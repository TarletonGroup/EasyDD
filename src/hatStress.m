function sigma = hatStress(uhat, nc, x, D, mx, my, mz, w, h, d, x0)

    % returns stress tensor at a point x=x0 by constructing B matrix B(x)
    % notes with mx=90 dx=6,dy=dz=1, agrees with S(1:5) agree with Abaqus to 1%
    % S(6) is approx zero so the error is huge!

    i = ceil(x0(1) / w);
    i = max(i, 1);
    j = ceil(x0(2) / h);
    j = max(j, 1);
    k = ceil(x0(3) / d);
    k = max(k, 1);
    p = i + (k - 1) * mx + (j - 1) * mx * mz;

    if any(x0 < 0) || p > mx * my * mz || isnan(p)
        % Sort of have to do this when using int_trapezium.m because real
        % segment will move out of domain during trial before it can be
        % remeshed?
        % fprintf('Node outside domain! See hatStress.m\n')
        sigma = zeros(3); % do something better like remeshing later...
        return
    end

    xc = zeros(3, 1);

    for i = 1:3
        xc(i) = 0.5 * (x(nc(p, 1), i) + x(nc(p, 7), i)); %xc is center of element p
    end

    a = w / 2; % element dx
    b = h / 2; % element dy
    c = d / 2; % element dz

    s1 = (x0(1) - xc(1)) / a;
    s2 = (x0(2) - xc(2)) / b;
    s3 = (x0(3) - xc(3)) / c;

    ds1dx = 1 / a;
    ds2dy = 1 / b;
    ds3dz = 1 / c;

    B = zeros(6, 24);
    dNds1 = zeros(8, 1);
    dNds2 = zeros(8, 1);
    dNds3 = zeros(8, 1);

    % 4. ----- .3
    %  |\       |\
    %  | \      | \
    % 1. -\---- .2 \
    %   \  \     \  \
    %    \ 8. ----\- .7
    %     \ |      \ |
    %      \|       \|
    %      5. ----- .6
    %  s2   s3
    %    \  ^
    %     \ |
    %      \|
    %       . -----> s1
    % redefine local corordinate system (s1,s2,s3) to have same orientation as
    % the global system (x,y,z) this should save calcualting Jij and inv(J). These shape function derivatives are built using the above geometry. They differ from the ones in `finiteElement3D.m` in that the ones found there are for a different FEM construction, i.e. building the mesh in the xy-plane and moving up the z-axis. The mesh builder in `finiteElement3D.m` builds the mesh in xz-planes travelling up the y-axis and so needs to define a Jacobian to transform those coordinates into those of the shape functions found there. Since these shape functions are custom for the geometry used, we don't need a Jacobian, we can simply use them to build B.
    pm1 = [-1 1 1 -1 -1 1 1 -1];
    pm2 = [1 1 1 1 -1 -1 -1 -1];
    pm3 = [-1 -1 1 1 -1 -1 1 1];
    %eg shape function a: Na = 1/8*(1+pm1(a)*s1)(1+pm2(a)*s2)(1+pm3(a)*s3)
    % dNa/ds1 = 1/8* pm1(a)*(1+pm2(a)*s2)(1+pm3(a)*s3)
    % dNa/dx = (dNa/ds1)(ds1/dx) where ds1/dx = 1/a
    % dN1/dx = 1/a(dNa/ds1) = -b*c*(1+s2)(1-s3)/(abc)

    for a = 1:8
        dNds1(a) = 1/8 * pm1(a) * (1 + pm2(a) * s2) * (1 + pm3(a) * s3);
        dNds2(a) = 1/8 * (1 + pm1(a) * s1) * pm2(a) * (1 + pm3(a) * s3);
        dNds3(a) = 1/8 * (1 + pm1(a) * s1) * (1 + pm2(a) * s2) * pm3(a);

        B(1, 3 * (a - 1) + 1) = dNds1(a) * ds1dx;
        B(2, 3 * (a - 1) + 2) = dNds2(a) * ds2dy;
        B(3, 3 * (a - 1) + 3) = dNds3(a) * ds3dz;

        B(4, 3 * (a - 1) + 1) = B(2, 3 * (a - 1) + 2);
        B(4, 3 * (a - 1) + 2) = B(1, 3 * (a - 1) + 1);

        B(5, 3 * (a - 1) + 1) = B(3, 3 * a);
        B(5, 3 * a) = B(1, 3 * (a - 1) + 1);

        B(6, 3 * (a - 1) + 2) = B(3, 3 * a);
        B(6, 3 * a) = B(2, 3 * (a - 1) + 2);
    end

    %
    %----------------------------------------------------------------------
    U = zeros(24, 1);

    for a = 1:8
        U(3 * a - 2) = uhat(3 * nc(p, a) - 2);
        U(3 * a - 1) = uhat(3 * nc(p, a) - 1);
        U(3 * a) = uhat(3 * nc(p, a));
    end

    sigmaA = D * (B * U); %

    sigma = zeros(3);
    sigma(1, 1) = sigmaA(1); %11
    sigma(2, 2) = sigmaA(2); %22
    sigma(3, 3) = sigmaA(3); % 33
    sigma(1, 2) = sigmaA(4); % 12
    sigma(1, 3) = sigmaA(5); % 13
    sigma(2, 3) = sigmaA(6); % 23
    sigma(2, 1) = sigma(1, 2);
    sigma(3, 1) = sigma(1, 3);
    sigma(3, 2) = sigma(2, 3);

end
