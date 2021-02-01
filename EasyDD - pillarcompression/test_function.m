%test NodalSurfForceLinearRectangle2.m and its mex version

for i=1:100
    
    x1=rand(1,3);
    x2=rand(1,3);
    x3=rand(1,3);
    x4=rand(1,3);
    x5=rand(1,3);
    x6=rand(1,3);
    b=rand(1,3);
    nu=rand(1,1);
    mu=rand(1,1);
    a=rand(1,1);
    t=rand(1,3);
    
    disp('Running Matlab version (Queyreau)...');
    tic;
    [fx3,fx4,fx5,fx6,ftot] = NodalSurfForceLinearRectangle2(x1,x2,x3,x4,x5,x6,b,mu,nu,a,t);
    time1=toc;
    
    disp('Running C version (Ferroni)...');
    tic;
    [fx3mex,fx4mex,fx5mex,fx6mex,ftotmex] = NodalSurfForceLinearRectangleMex(x1,x2,x3,x4,x5,x6,b,mu,nu,a,t);
    time2=toc;
   
    speedup = time1/time2;
    disp('Speed-up of MEX function');
    disp(speedup);
    
    isequalAbs = @(x,y,tol) ( abs(x-y) <= tol );
    
    test1 = all(isequalAbs(fx3,fx3mex',10e-14));
    test2 = all(isequalAbs(fx4,fx4mex',10e-14));
    test3 = all(isequalAbs(fx5,fx5mex',10e-14));
    test4 = all(isequalAbs(fx6,fx6mex',10e-14));
    
    if (test1 && test2 && test3 && test4 == 1)
        disp('Results match! :) ');
    else
        disp('Results dont match! :( ');
    end
end
    
    