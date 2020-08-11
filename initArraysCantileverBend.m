function [Fend, U_bar, t, simTime, dt, curstep] = ...
    initArraysCantileverBend(dt0, maxSteps)
%=========================================================================%
% Initialises simulation variables.
%
% Daniel Celis Garza, Aug 2020
%-------------------------------------------------------------------------%
% Inputs
% dt0 := initial time step.
% maxSteps := maximum number of steps in a simulation.
%-------------------------------------------------------------------------%
% Outputs
% Fend := array of forces on FE nodes at ever time step
% U_bar := displacements on FE nodes
% t := array of times at each step
% simTime : = current simulation time
% dt := current dt
% curstep := current step
%=========================================================================%

Fend = zeros(maxSteps*3,1); 
U_bar = zeros(maxSteps*3,1);
t = zeros(maxSteps*3,1); 
simTime = 0;
dt = dt0;
curstep = 0;

end