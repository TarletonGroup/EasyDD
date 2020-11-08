function [curstep, simTime] = updateTime(curstep, simTime, dt)
    curstep = curstep + 1;
    simTime = simTime + dt;
end
