function [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = collideNodesAndSegments(...
    rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, ...
    u_hat, curstep, ...
    matpara, mods, flags, FEM, Bcoeff)
    %% Extraction
    
    % matpara:
    rann = matpara.rann;

    % flags:
    docollision = flags.docollision;

    %% Function

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
                    [rnnew, linksnew, ~, ~, fsegnew, colliding_segments] = collision(...
                        rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, ...
                        u_hat, floop, n1s1, n2s1, n1s2, n2s2, s1, s2, segpair, ...
                        matpara, mods, flags, FEM, Bcoeff);

                    %removing links with effective zero Burgers vectors
                    [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = cleanupsegments(rnnew, linksnew, fsegnew);
                end

            end

        end

    end

end
