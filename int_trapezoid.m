function [rn,vn,dt,fn,fseg]=int_trapezoid(rn,dt,dt0,MU,NU,a,Ec,links,connectivity,...
    rmax,rntol,mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d)
%implicit numerical integrator using the Trapezoid method
%dt: suggested timestep (usually from previous iteration)
%dt0: maximum allowed timestep

%dummy variable
t=0;
rnold=rn;

%scramble rn into a single column vector
rnvec=[rn(:,1);rn(:,2);rn(:,3)]; flag=rn(:,4);

%Backward Euler
rnvec0=rnvec;

[vnvec0,fn,fseg]=drndt(rnvec0,flag,MU,NU,a,Ec,links,connectivity,...
        mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);
%dt=1/max(1/dt0,max(vnvec0)/rmax);
%dt=dt0;

dt1=dt;
maxiter=1;
convergent=0;

% This algorithm [Cai & Bulatov, Algorithm 10.2, pg. 216] is a simple
% routine to dynamically change the time-step. Whenever the difference
% between the predicted and corrected positions for any node exceeds
% threshold err, the timestep is decreased by half and new positions for
% all nodes are recomputed with the reduced time step.
%   1. Initialize time step dt = dt_max
%   2. dt_0 = dt
%   3. Compute predictor rn_P(t+dt) and corrector rn(t+dt) from forward
%   Euler and trapezoid method respectively.
%   4. If max(||rn_P(t+dt) -  rn(t+dt)||) > e, reduce timestep by half
%   5. t = t + dt
%   6. If dt = dt_0, increase the step to dt=min(1.2dt, dt_max)
%   7. Return to 2, unless total number of cycles is reached

while(~convergent)
    
    rnvec1=rnvec0+vnvec0*dt; %Euler forward method [Cai & Bulatov, eq. 10.43]
    if isempty(rnvec1)
        convergent=1;
    end
% ET - rnvec1 can contain a node outside the domain! 
    for iter=1:maxiter,
        [vnvec,fn,fseg]=drndt(rnvec1,flag,MU,NU,a,Ec,links,connectivity,...
            mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);
        %err=rnvec1-rnvec0-vnvec.*dt;          %backward Euler
        err=rnvec1-rnvec0-(vnvec+vnvec0)/2*dt; %trapzoid
        errmag=max(abs(err));
        %disp(sprintf('iter=%d err=%e',iter,errmag));
        if(errmag<rntol)
            convergent=1;
            break;
        else
            rnvec1=rnvec1-err;
        end
    end
    
    if(convergent)
        break;
    else
        dt=dt/2;
    end
end

%unscramble rn and vn vectors
rn=[reshape(rnvec1,length(rnvec1)/3,3),flag];
vn=reshape(vnvec,length(vnvec)/3,3); % trapezoidal rule modification

%When no time step reduction is necessary, the algorithm attempts to
%increase dt by 20% for the next cycle, but not to exceed a preset value of
%dtmax. 

%automatically adjust time step
if isempty(errmag)
    dt=dt0;
elseif((dt==dt1)&&(iter==1))
    maxchange=1.2;
    exponent=20;
    factor=maxchange*(1/(1+(maxchange^exponent-1)*(errmag/rntol)))^(1/exponent);
    dt=min(dt1*factor,dt0);
end