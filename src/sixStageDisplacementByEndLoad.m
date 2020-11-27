function [sign_u_dot, u_dot, sign_f_dot, f_dot, u_tilda_0, u, u_hat, u_tilda] = sixStageDisplacementByEndLoad(...
        sign_u_dot, u_dot, sign_f_dot, f_dot, u_tilda_0, u, f_bar, f_hat, f_tilda, u_bar, ...
        u_hat, u_tilda, r_hat, simTime, totalSimTime, curstep, loadingFunctionArgStruct)

    u_dot_0 = loadingFunctionArgStruct.u_dot_0;
    u_bar_crit = loadingFunctionArgStruct.u_bar_crit;
    scaleFactor = loadingFunctionArgStruct.scaleFactor;

    if u_bar < u_bar_crit(1)
        u_dot = u_dot_0;
    elseif u_bar_crit(1) < u_bar && u_bar < u_bar_crit(2)
        u_dot = u_dot_0 * scaleFactor(1);
    elseif u_bar_crit(2) < u_bar && u_bar < u_bar_crit(3)
        u_dot = u_dot_0 * scaleFactor(2);
    elseif u_bar_crit(3) < u_bar && u_bar < u_bar_crit(4)
        u_dot = u_dot_0 * scaleFactor(3);
    elseif u_bar_crit(4) < u_bar && u_bar < u_bar_crit(5)
        u_dot = u_dot_0 * scaleFactor(4);
    else
        u_dot = u_dot_0 * scaleFactor(5);
    end

end
