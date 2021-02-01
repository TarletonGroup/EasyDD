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
    
    fprintf('Running Matlab version (Queyreau)...\n');
    tic;
    [fx3,fx4,fx5,fx6,ftot] = NodalSurfForceLinearRectangle2(x1,x2,x3,x4,x5,x6,b,mu,nu,a);
    time1=toc;
    
    fprintf('Running C version (Ferroni)...\n');
    tic;
    [fx3mex,fx4mex,fx5mex,fx6mex,ftotmex] = NodalSurfForceLinearRectangleMex(x1,x2,x3,x4,x5,x6,b,mu,nu,a);
    time2=toc;
  
    speedup = time1/time2;
    fprintf('Speed-up of manually-written MEX function\n');
    fprintf('%d\n',speedup);
    
     tic;
     fprintf('Running C version (Matlab CODER 2.7)...\n');
     [fx3mex2,fx4mex2,fx5mex2,fx6mex2,ftotmex2] = NodalSurfForceLinearRectangleMex(x1,x2,x3,x4,x5,x6,b,mu,nu,a);
     time3=toc;
    
    fprintf('Speed-up of automatic MEX function\n');
    fprintf('%d\n',time1/time3);
    isequalAbs = @(x,y,tol) ( abs(x-y) <= tol );
    
%     test1 = all(isequalAbs(fx3,fx3mex',10e-14));
%     test2 = all(isequalAbs(fx4,fx4mex',10e-14));
%     test3 = all(isequalAbs(fx5,fx5mex',10e-14));
%     test4 = all(isequalAbs(fx6,fx6mex',10e-14));
    
    test1b = all(isequalAbs(fx3,fx3mex2',10e-14));
    test2b = all(isequalAbs(fx4,fx4mex2',10e-14));
    test3b = all(isequalAbs(fx5,fx5mex2',10e-14));
    test4b = all(isequalAbs(fx6,fx6mex2',10e-14));
   
%     if (test1 && test2 && test3 && test4 == 1)
%         fprintf('Results between Matlab and C (Ferroni) match! :)\n\n');
%     else
%         fprintf('Results dont match! :(\n\n');
%     end
     
    if (test1b && test2b && test3b && test4b == 1)
        fprintf('Results between Matlab and C (Coder 2.7) match! :)\n\n');
    else
        fprintf('Results dont match! :(\n\n');
    end

end