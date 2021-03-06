function [rn, vn, dt, fn, fseg] = int_trapezoid(rn, dt, dt0, dtMin, MU, NU, a, Ec, links, connectivity, ...
        rmax, rntol, mobility, vertices, rotMatrix, uhat, nc, xnodes, D, mx, my, mz, w, h, d, Bcoeff, CUDA_flag)

    %Implicit numerical integrator using the Euler-trapezoid method adapted
    %from [Cai & Bulatov, Algorithm 10.2, pg. 216]. The tiemstep size is
    %controlled so that it can increase suitably quickly and not decrease
    %excessively whilst remaining within acceptable tolerence limits
    %Written by B.Bromage and D.Celis-Garza 05/11/2018

    %Convert rn into a single column of coordinates and store the node flags
    rnvec0 = [rn(:, 1); rn(:, 2); rn(:, 3)]; flag = rn(:, 4);

    %Calculate the current nodal velocities
    [vnvec0, fn, fseg] = drndt(rnvec0, flag, MU, NU, a, Ec, links, connectivity, ...
        mobility, vertices, rotMatrix, uhat, nc, xnodes, D, mx, my, mz, w, h, d, Bcoeff, CUDA_flag);

    maxiter = 10; %Maximum number of times the timestep can be increased
    counter = 1; %Counter variable for the while loop
    dt_old = 0; %Variable for the maximumum acceptable timestep
    maxchange = 1.2; %Maximum factor the timestep can be increased by
    exponent = 20; %Variable that controls the degree that the calculated error affects the cange in timestep
    dt_old_good = 0; %Logical which flags whether an acceptable timestep has been calculated
    convergent = 0; %Logical for the operation of the while loop
    while (~convergent)
        rnvec1 = rnvec0 + vnvec0 * dt; %Euler forward method [Cai & Bulatov, eq. 10.43]

        if isempty(rnvec1)%If there are no dislocations then use the maximum possible timestep
            convergent = 1;
        end

        %Calculate the nodal velocities for the next timestep accordin to
        %Euler forward method
        [vnvec1, fn, fseg] = drndt(rnvec1, flag, MU, NU, a, Ec, links, connectivity, ...
            mobility, vertices, rotMatrix, uhat, nc, xnodes, D, mx, my, mz, w, h, d, Bcoeff, CUDA_flag);
        distvec = rnvec1 - rnvec0;
        
        distmag = max(sum(reshape(distvec, [], 3) .* reshape(distvec, [], 3), 2));%max(abs(distvec)); %Largest distance any node moves
%         distmag = quantile(abs(distvec), 0.99);
        err = distvec - (vnvec1 + vnvec0) / 2 * dt; %Euler-trapezoid method
        errmag = max(sum(reshape(err, [], 3) .* reshape(err, [], 3), 2));%max(abs(err)); %Largest difference (error) calculated by the Euler-trapezoid method
%         errmag = quantile(abs(err), 0.99);
        if isempty(errmag)%If the Euler-trapzoid method yields no results use maximum time step and end loop
            dt = dt0;
            break
        end

        if (errmag < rntol) && (distmag < rmax)%If error and max distance move are in acceptable limits
            dt_old = dt; %Store current timestep as maximum acceptable timestep
            factor = maxchange * (1 / (1 + (maxchange^exponent - 1) * (errmag / rntol)))^(1 / exponent);
            dt = min(dt * factor, dt0); %Increase timestep depending on the magnitude of the error
            dt_old_good = 1; %Flag acceptable timestep calculated
            counter = counter + 1; %Proceed to next iteration
        elseif dt_old_good == 1%If current timestep is too large, use largest acceptable timestep
            dt = dt_old;
            counter = maxiter;
        elseif dt < dtMin
            counter = counter + 1;
            dt = dtMin;
        elseif dt == dtMin
            counter = maxiter + 1;
        else
            dt = dt / 2; %If no acceptable timestep has been calculated, halve timestep and try again
        end

        if counter > maxiter || dt == dt0%End loop if maximum number of iterations is reached or curren timestep is maximum timestep
            convergent = 1;
        end

    end

    %Rearrange rnvec and vnvec into rn and vn matrices
    rn = [reshape(rnvec1, length(rnvec1) / 3, 3), flag];
    vn = reshape(vnvec1, length(vnvec1) / 3, 3);
end
