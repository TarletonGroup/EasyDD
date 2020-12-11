function [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = collideNodesAndSegments(docollision, ...
        rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, rann, MU, NU, a, Ec, mobility, vertices, rotMatrix, ...
        u_hat, nc, xnodes, D, mx, mz, w, h, d, lmin, CUDA_flag, Bcoeff, curstep)

    if (docollision)

        % Collision detection and handling
        colliding_segments = 1;

        while colliding_segments == 1
            [colliding_segments, n1s1, n2s1, n1s2, n2s2, floop, s1, s2, segpair] = CollisionCheckerMex(rnnew(:, 1), rnnew(:, 2), rnnew(:, 3), rnnew(:, end), ...
                rnnew(:, 4), rnnew(:, 5), rnnew(:, 6), linksnew(:, 1), linksnew(:, 2), connectivitynew, rann);

            if colliding_segments == 1%scan and update dislocation structure.

                if floop == 1
                    fprintf("Step %d. Unconnected links found. Links %d and %d are colliding.\n", curstep, s1, s2)
                elseif floop == 2
                    fprintf("Step %d. Links %d and %d colliding by hinge condition.\n", curstep, s1, s2)
                end

                if colliding_segments == 1
                    [rnnew, linksnew, ~, ~, fsegnew, colliding_segments] = collision(rnnew, linksnew, connectivitynew, ...
                        linksinconnectnew, fsegnew, rann, MU, NU, a, Ec, mobility, vertices, rotMatrix, u_hat, nc, xnodes, ...
                        D, mx, mz, w, h, d, floop, n1s1, n2s1, n1s2, n2s2, s1, s2, segpair, lmin, CUDA_flag, Bcoeff);

                    %removing links with effective zero Burgers vectors
                    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = cleanupsegments(rnnew, linksnew, fsegnew);
                end

            end

        end

    end

end
