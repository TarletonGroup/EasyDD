function [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = collideAndSeparateNodesAndSegments(docollision, doseparation, ...
        rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, rann, MU, NU, a, Ec, mobility, vertices, rotMatrix, ...
        u_hat, nc, xnodes, D, mx, mz, w, h, d, lmin, CUDA_flag, Bcoeff, curstep)

    if (docollision)

        % Collision detection and handling
        colliding_segments = 1;

        while colliding_segments == 1
            [colliding_segments, n1s1, n2s1, n1s2, n2s2, floop, s1, s2, segpair] = CollisionCheckerMex(rnnew(:, 1), rnnew(:, 2), rnnew(:, 3), rnnew(:, end), ...
                rnnew(:, 4), rnnew(:, 5), rnnew(:, 6), linksnew(:, 1), linksnew(:, 2), connectivitynew, rann,linksnew(:,3),linksnew(:,4),linksnew(:,5),linksnew(:,6),linksnew(:,7),linksnew(:,8));

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
                    
                    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = separation(doseparation, rnnew, ...
                        linksnew, connectivitynew, linksinconnectnew, fsegnew, mobility, rotMatrix, MU, NU, a, Ec, ...
                        2 * rann, vertices, u_hat, nc, xnodes, D, mx, mz, w, h, d, CUDA_flag, Bcoeff);

                end

            end

        end

    end

end
