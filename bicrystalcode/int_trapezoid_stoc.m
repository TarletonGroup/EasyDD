function [rn,vn,dt,fn,fseg]=int_trapezoid_stoc(rn,dt,dt0,MU,NU,a,Ec,links,connectivity,...
    ~,rntol,mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d)
% implicit numerical integrator using the Trapezoid method
% dt: suggested timestep (usually from previous iteration)
% dt0: maximum allowed timestep
%
% includes simplectic algorithm for stochastic noise in langevin simulation
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NB NOISE IS NOT CONSTRAINED BASED ON CRYSTALLOGRAPHY!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global LangevinDiffusionCoeff;

%scramble rn into a single column vector
rnvec=[rn(:,1);rn(:,2);rn(:,3)]; flag=rn(:,4);

%Generate thermal noise (nodal based), from normal distribution with mean=0
%scaled with sqrt(2*D), where D is the diffusion coefficient.

% S.L.Dudarev et al./Journal of Nuclear Materials 455 (2014) 16-20
% D(T) = (Dbar/2rho) exp(-Ea/kT), where Dbar = 5.7x10^13 nm^3/s, rho is the
% loop radius (in nm units), Ea=1.3eV is the activation for loop migration,
% and Eb approx. 0.4eV
%what is D in DDLab units? Scales with loop radius apparently?
vn_langevin = sqrt(2*LangevinDiffusionCoeff).*randn(size(rnvec,1),size(rnvec,2));

%Backward Euler
%rnvec0=rnvec; 
rnvec0 = rnvec + vn_langevin*0.5*dt; %Step 1 of simplectic algorithm

%predict velocities from elastic interactions (predictor step)
[vnvec0,fn,fseg]=drndt(rnvec0,flag,MU,NU,a,Ec,links,connectivity,...
        mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);

    
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
                             %This is also Step 2 of simplectic algorithm.
    
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
        %Step 3 of the simplectic algorithm, i.e. 2nd part of noise
        vnvec = vnvec + vn_langevin*0.5*dt;
        break;
    else
        dt=dt/2;
        %Update Step 1 from simplectic algorithm...
        %This is important since dt needs to be the same for all three
        %steps of the simplectic algorithm. So, if we incorporate this in
        %the Euler forward method, we need to update rnvec0 to account for
        %the new dt being applied to the noise part.
        rnvec0 = rnvec + vn_langevin*0.5*dt;
        [vnvec0,fn,fseg]=drndt(rnvec0,flag,MU,NU,a,Ec,links,connectivity,...
        mobility,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d);
    end
end

%unscramble rn and vn vectors
rn=[reshape(rnvec1,length(rnvec1)/3,3),flag];
vn=reshape(vnvec,length(vnvec)/3,3); % trapezoidal rule modification

%When no time step reduction is necessary, the algorithm attempts to
%increase dt by 20% for the next cycle, but not to exceed a preset value of
%dtmax. 

%automatically adjust time step
if((dt==dt1)&&(iter==1))
    maxchange=1.2;
    exponent=20;
    factor=maxchange*(1/(1+(maxchange^exponent-1)*(errmag/rntol)))^(1/exponent);
    dt=min(dt1*factor,dt0);
end