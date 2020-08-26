function [Fsim, Usim, t] = updateForceDispTime(Fsim, Usim, t, f_out, u_out, simTime, curstep, simType)

    if simType == 1
        Fsim(curstep + 1) = f_out;
        Usim(curstep + 1) = u_out;
        t(curstep + 1) = simTime;
    end

end
