function [output]=plot_int_err(rn,dt0,MU,NU,a,Ec,links,connectivity,mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d)

rnvec0=[rn(:,1);rn(:,2);rn(:,3)]; flag=rn(:,4);

[vnvec0,~,~]=drndt(rnvec0,flag,MU,NU,a,Ec,links,connectivity,...
        mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);
    
    errmag=zeros([1001 1]);
    dtmag=zeros([1001 1]);
    
    for i=1:1001
        dt=(i-1)*dt0/(1000);
        dtmag(i)=dt;
        rnvec1=rnvec0+vnvec0*dt;
        [vnvec,~,~]=drndt(rnvec1,flag,MU,NU,a,Ec,links,connectivity,...
            mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);
        err=rnvec1-rnvec0-(vnvec+vnvec0)/2*dt;
        errmag(i)=max(abs(err));
    end
    
    plot(dtmag,errmag)
    
    plotHandle=gcf;
    
    output=plotHandle;
    
    