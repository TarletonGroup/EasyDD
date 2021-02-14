function [rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew] = collideAndSeparateNodesAndSegments(docollision, doseparation, ...
        rnnew, linksnew, connectivitynew, linksinconnectnew, fsegnew, rann, MU, NU, a, Ec, mobility, vertices, rotMatrix, ...
        u_hat, nc, xnodes, D, mx, my, mz, w, h, d, lmin, CUDA_segseg_flag, Bcoeff, curstep)

    if (docollision)

        % Collision detection and handling
        colliding_segments = 1;
        
        s1Skip = [];
        s2Skip = [];
        while colliding_segments == 1
            [colliding_segments, n1s1, n2s1, n1s2, n2s2, floop, s1, s2, segpair] = CollisionCheckerMex(rnnew(:, 1), rnnew(:, 2), rnnew(:, 3), rnnew(:, end), ...
                rnnew(:, 4), rnnew(:, 5), rnnew(:, 6), linksnew(:, 1), linksnew(:, 2), connectivitynew, rann,linksnew(:,3),linksnew(:,4),linksnew(:,5),linksnew(:,6),linksnew(:,7),linksnew(:,8), s1Skip, s2Skip);

            if colliding_segments == 1%scan and update dislocation structure.

                if floop == 1
                    fprintf("Step %d. Unconnected links found. Links %d and %d are colliding.\n", curstep, s1, s2)
                elseif floop == 2
                    fprintf("Step %d. Links %d and %d colliding by hinge condition.\n", curstep, s1, s2)
                end
                
                numNode = size(rnnew, 1);
                numSeg = size(linksnew, 1);

                if colliding_segments == 1
                    [rnnewTmp, linksnewTmp, ~, ~, fsegnewTmp, colliding_segments] = collision(rnnew, linksnew, connectivitynew, ...
                        linksinconnectnew, fsegnew, rann, MU, NU, a, Ec, mobility, vertices, rotMatrix, u_hat, nc, xnodes, ...
                        D, mx, my, mz, w, h, d, floop, n1s1, n2s1, n1s2, n2s2, s1, s2, segpair, lmin, CUDA_segseg_flag, Bcoeff);

                    %removing links with effective zero Burgers vectors
                    [rnnewTmp, linksnewTmp, connectivitynewTmp, linksinconnectnewTmp, fsegnewTmp] = cleanupsegments(rnnewTmp, linksnewTmp, fsegnewTmp);
                    
                    [rnnewTmp, linksnewTmp, connectivitynewTmp, linksinconnectnewTmp, fsegnewTmp] = separation(doseparation, rnnewTmp, ...
                        linksnewTmp, connectivitynewTmp, linksinconnectnewTmp, fsegnewTmp, mobility, rotMatrix, MU, NU, a, Ec, ...
                        2 * rann, vertices, u_hat, nc, xnodes, D, mx, my, mz, w, h, d, CUDA_segseg_flag, Bcoeff); 
                    numNodeTmp = size(rnnewTmp, 1);
                    numSegTmp = size(linksnewTmp, 1);
                    if numNodeTmp == numNode && numSegTmp == numSeg % || numNodeTmp == numNode - 1 && numSegTmp == numSeg - 1
                        s1Skip = [s1Skip; s1];
                        s2Skip = [s2Skip; s2];
                        
                        if size(s1Skip, 1) > 1 && (s1 == s1Skip(end) && s2 == s2Skip(end) || s2 == s1Skip(end) && s1 == s2Skip(end))
                            rnnew = rnnewTmp;
                            linksnew = linksnewTmp;
                            connectivitynew = connectivitynewTmp;
                            linksinconnectnew = linksinconnectnewTmp;
                            fsegnew = fsegnewTmp;
                            s1Skip = [];
                            s2Skip = [];
                        end
                        
                        continue
                    else
                        rnnew = rnnewTmp;
                        linksnew = linksnewTmp;
                        connectivitynew = connectivitynewTmp;
                        linksinconnectnew = linksinconnectnewTmp;
                        fsegnew = fsegnewTmp;
                        s1Skip = [];
                        s2Skip = [];
                    end

                end

            end

        end

    end

end
