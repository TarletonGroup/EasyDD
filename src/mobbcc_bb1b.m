function [vn, fn] = mobbcc_bb1b(fseg, rn, links, connectivity, nodelist, conlist, Bcoeff, rotMatrix)
    %mobility law function (model: BCC0)

    Bscrew = Bcoeff.screw;
    Bedge = Bcoeff.edge;
    Beclimb = Bcoeff.climb;
    rotateCoords = false;

    if ~isempty(rotMatrix)
        rotateCoords = true;
        rn(:, 1:3) = rn(:, 1:3) * rotMatrix;
        fseg(:, 1:3) = fseg(:, 1:3) * rotMatrix;
        fseg(:, 4:6) = fseg(:, 4:6) * rotMatrix;
        links(:, 3:5) = links(:, 3:5) * rotMatrix;
        links(:, 6:8) = links(:, 6:8) * rotMatrix;
    end

    %numerical tolerance
    tol = 1e-7;
    %list of permitted slip planes
    planelist = (1 / sqrt(2)) * [1 1 0;
                            1 0 1;
                            0 1 1;
                            1 -1 0;
                            1 0 -1;
                            0 1 -1];

    % length of the nodelist for which the velocity will be calculated
    L1 = size(nodelist, 1);
    % if no nodelist is given then the nodelist becomes the whole node population
    % this portion of the code sets up that nodelist along with the connlist
    % that contains all of the nodal connections
    if L1 == 0
        L1 = size(rn, 1);
        nodelist = linspace(1, L1, L1)';
        [L2, L3] = size(connectivity);
        conlist = zeros(L2, (L3 - 1) / 2 + 1);
        conlist(:, 1) = connectivity(:, 1);

        for i = 1:L2
            connumb = conlist(i, 1);
            conlist(i, 2:connumb + 1) = linspace(1, connumb, connumb);
        end

    end

    % now cycle through all of the nodes for which the velocity must be calculated

    vn = zeros(L1, 3);
    fn = zeros(L1, 3);

    for n = 1:L1
        n0 = nodelist(n); %n0 is the nodeid of the nth node in nodelist
        numNbrs = conlist(n, 1); %numNbrs is the number of connections for node n0 in conlist
        Btotal = zeros(3, 3);

        for i = 1:numNbrs
            ii = conlist(n, i + 1); % connectionid for this connection
            linkid = connectivity(n0, 2 * ii);
            posinlink = connectivity(n0, 2 * ii + 1);
            n1 = links(linkid, 3 - posinlink);

            if rn(n1, end) == 67
                continue
            end

            rt = rn(n1, 1:3) - rn(n0, 1:3); % calculate the length of the link and its line direction
            L = norm(rt);

            if L > 0.0% only calculate drag for non zero length links
                fsegn0 = fseg(linkid, 3 * (posinlink - 1) + (1:3));
                %             if norm(fsegn0)/L<0.0015  % for implementing P-N stress
                %                 fsegn0=[0,0,0];
                %             end
                fn(n, :) = fn(n, :) + fsegn0; % nodeid for the node that n0 is connected to
                burgv = links(connectivity(n0, 2 * ii), 3:5); % burgers vector of the link
                mag = norm(burgv);
                checkv = abs(burgv);
                checkv = checkv / min(checkv);
                check111 = sum(abs(1 - checkv(:)) < tol); % check if burgers vector is <111> type
                linedir = rt ./ L; % normalise line direction
                checkmag = abs(sqrt(3) * 0.5 - mag); % check if burgers vector is single or compound/super dislocation

                if check111 == 3 && checkmag <= eps% if burgeers vector is a single <111> type
                    costh2 = (linedir * burgv')^2 / (burgv * burgv'); % (lhat.bhat)^2 = cos^2(theta)% calculate how close to screw the link is
                    sinth2 = 1 - costh2;

                    if sinth2 > tol% not pure screw segment
                        dotprods = planelist * burgv' / mag; % find slip planes burgers vector can exist in
                        slipplanes = planelist(abs(dotprods) < eps, :);
                        dotprods2 = slipplanes * linedir'; % find available slipplane closest to dislocation habit plane
                        cosdev = dotprods2(abs(dotprods2) == min(abs(dotprods2)));

                        if size(cosdev, 1) > 1% if multiple slip planes are listed, choose one. This makes little difference since mobility will be low if this is the case
                            cosdev = cosdev(1);
                        end

                        ndir = slipplanes(dotprods2 == cosdev, :);

                        if size(ndir, 1) > 1
                            ndir = ndir(1, :);
                        end

                        mdir = cross(ndir, linedir); % find glide direction
                        clinedir = cross(mdir, ndir); % find projection of line direction into the slip plane
                        linecos2 = (linedir * clinedir')^2; % find how out of plane the line direction is for interpolation of line drag
                        linesin2 = 1 - linecos2;
                        Bglide = 1 / sqrt((1 / Bedge^2) * sinth2 + (1 / Bscrew^2) * costh2); % Eqn (112) from Arsenlis et al 2007 MSMSE 15 553
                        Bglide = 1 / sqrt((1 / Beclimb^2) * linesin2 + (1 / Bglide^2) * linecos2); % Interpolate glide and line coefficients for how screw type and how out of plane they are
                        Bline2 = 1 / sqrt((1 / Bscrew^2) * sinth2 + (1 / Bedge^2) * costh2);
                        Bline2 = 1 / sqrt((1 / Beclimb^2) * linesin2 + (1 / Bline2^2) * linecos2);
                        Btotal = Btotal + mag .* ((0.5 * L) .* ((Bglide) .* (mdir' * mdir) + (Beclimb) .* (ndir' * ndir) + (Bline2) .* (clinedir' * clinedir))); % construct drag tensor for this segment
                    else % pure screw segment

                        if norm(fsegn0) > eps% find normalised direction of force on the segment for node in question
                            fnorm = fsegn0 / norm(fsegn0);
                        else
                            fnorm = [0 0 0];
                        end

                        dotprods = planelist * burgv' / mag; % find slip planes burgers vector exists in
                        slipplanes = planelist(abs(dotprods) == min(abs(dotprods)), :);
                        dotprods2 = slipplanes * fnorm'; % find closest slip plane that force vector lies in
                        cosdev = dotprods2(abs(dotprods2) == min(abs(dotprods2)));

                        if size(cosdev, 1) > 1
                            cosdev = cosdev(1);
                        end

                        ndir = slipplanes(dotprods2 == cosdev, :); % assign slip plane normal
                        fsegn1 = fseg(linkid, 3 * (2 - posinlink) + (1:3)); % repeat process for force component acting on other node

                        if norm(fsegn1) > eps
                            fnorm_alt = fsegn1 / norm(fsegn1);
                        else
                            fnorm_alt = [0 0 0];
                        end

                        dotprods2_alt = slipplanes * fnorm_alt';
                        cosdev_alt = dotprods2_alt(abs(dotprods2_alt) == min(abs(dotprods2_alt)));

                        if size(cosdev_alt, 1) > 1
                            cosdev_alt = cosdev_alt(1);
                        end

                        if size(cosdev_alt, 1) > 1
                            cosdev_alt = cosdev_alt(1);
                        end

                        if ~isequal(size(cosdev_alt), [1 1]) ||~isequal(size(slipplanes, 1), size(dotprods2_alt, 1)) ||~isequal(size(dotprods2_alt, 2), 1)% check calculations are performed correctly
                            disp('YDFUS, see line 151 of mobbcc_bb1b')
                        end

                        ndir_alt = slipplanes(dotprods2_alt == cosdev_alt, :);

                        if size(ndir, 1) == 1 && size(ndir_alt, 1) == 1% check if both force components are acting in the same plane
                            planecheck = 1 - (ndir * ndir_alt');
                        end

                        if norm(fsegn0) < eps &&~norm(fsegn1) < eps% if one force component is zero, use the other to define the slip plane
                            ndir = ndir_alt;
                        elseif norm(fsegn1) < eps &&~norm(fsegn0) < eps
                            ndir_alt = ndir;
                        end

                        if size(ndir, 1) > 1 || size(ndir_alt, 1) > 1 || planecheck > eps% if the slip plane is uncertain
                            ndir = ndir(1, :); % choose an arbitrary slip plane for construction of drag matrix
                            mdir = cross(ndir, linedir);

                            if abs(cosdev) < eps
                                Btotal = Btotal + mag .* ((0.5 * L) .* ((Bscrew) .* (mdir' * mdir) + (Bscrew) .* (ndir' * ndir) + (Bedge) .* (linedir' * linedir))); % allow free glide in all directions if force is negligible, allowing other segments to dictate behaviour
                            else
                                Btotal = Btotal + mag .* ((0.5 * L) .* ((Beclimb) .* (mdir' * mdir) + (Beclimb) .* (ndir' * ndir) + (Bedge) .* (linedir' * linedir))); % restrict glide if force is significant, preventing undesired cross slip
                            end

                        else
                            mdir = cross(ndir, linedir); % find glide direction
                            cosdev2 = cosdev * cosdev; % find how out of plane the applied force is
                            cosratio = max(0, 1 - 4 * cosdev2); % calibrate interpolation so that glide is easiest on recognised glide planes and most difficult halfway between, cosratio must not be negative
                            sinratio = 1 - cosratio;
                            Bglide = 1 / sqrt((1 / Beclimb^2) * sinratio + (1 / Bscrew^2) * cosratio); % interpolate both glide and climb drag coefficients to restrict motion as force becomes more out of plane, while enabling a smooth transition of behaviour between slip systems
                            Bsclimb = 1 / sqrt((1 / Beclimb^2) * cosratio + (1 / Bscrew^2) * sinratio);
                            Btotal = Btotal + mag .* ((0.5 * L) .* ((Bglide) .* (mdir' * mdir) + (Bsclimb) .* (ndir' * ndir) + (Bedge) .* (linedir' * linedir))); % construct drag tensor as normal
                        end

                    end

                else %if this is not a viable slip system
                    costh2 = (linedir * burgv')^2 / (burgv * burgv'); % (lhat.bhat)^2 = cos^2(theta)% calculate how close to screw the link is
                    sinth2 = 1 - costh2;
                    Bline2 = 1 / sqrt((1 / Bscrew^2) * sinth2 + (1 / Bedge^2) * costh2); % Adapted from Eqn (112) from Arsenlis et al 2007 MSMSE 15 553
                    Btotal = Btotal + mag .* ((0.5 * L) .* ((Beclimb) .* eye(3) + (Bline2 - Beclimb) .* (linedir' * linedir))); % restrict motion in all directions except the line direction which is determined by how screw type te segment is
                end

            end

        end

        if norm(Btotal) < eps%if there is no drag on the node, make its velocity zero
            vn(n, :) = [0 0 0];
        elseif rcond(Btotal) < 1e-15%if the drag tensor is poorly conditioned use special inversion protocol
            Btotal_temp = Btotal + 1e-6 * max(max(abs(Btotal))) * eye(3); % perturb drag tensor
            Btotal_temp2 = Btotal - 1e-6 * max(max(abs(Btotal))) * eye(3);
            vn_temp = (Btotal_temp \ fn(n, :)')'; % estimate velocity using perturbations
            vn_temp2 = (Btotal_temp2 \ fn(n, :)')';
            vn(n, :) = 0.5 * (vn_temp + vn_temp2); % use mean of estimated velocities
        else
            vn(n, :) = (Btotal \ fn(n, :)')'; % Btotal was well conditioned so just take the inverse
        end

        if any(any(isnan(vn))) || any(any(~isreal(vn)))% ensure no non-physical velocities exist
            disp('YDFUS, see line 212 of mobbcc_bb1b')
        end

    end
    
    if rotateCoords
        vn = vn * rotMatrix';
        fn = fn * rotMatrix';
    end
end
