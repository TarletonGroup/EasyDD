function [rn,links] = cross_slip_BCCrotate5(fseg,rn,links,connectivity,MU,NU,a,Ec,vertices,uhat,nc,xnodes,D,mx,mz,w,h,d,areamin,L_cross_crit)

%[rn,links] = cross_slip_BCCrotate5(fseg,rn,links,connectivity,[],[],curstep,MU,NU,a,Ec,0,appliedstress);

thetacrit = 5/ 180.0 * 3.14159;

reset_pos_n1 = 0;
reset_pos_n2 = 0;
reset_pos_n = 0;
reset_plane_seg1 = 0;
reset_plane_seg2 = 0;

rotationBCC=eye(1);

rncrystal(:,1:3) = rn(:,1:3)*rotationBCC';
rncrystal(:,4) = rn(:,4);
fsegcrystal(:,1:3) = fseg(:,1:3)*rotationBCC';
fsegcrystal(:,4:6) = fseg(:,4:6)*rotationBCC';
linkscrystal(:,3:5) = links(:,3:5)*rotationBCC';
linkscrystal(:,6:8) = links(:,6:8)*rotationBCC';
linkscrystal(:,1:2) = links(:,1:2);
linkscrystal(find(abs(linkscrystal)<1E-4)) = 0;

eps = 1E-10;


sthetacrit = sin(thetacrit);
s2thetacrit = sthetacrit * sthetacrit;
MU = 1;

noiseFactor=1e-03;
weightFactor=1.0;

[longscrewsegs,longscrewnodes] = findscrew(rncrystal,linkscrystal,connectivity,s2thetacrit);

nodelist = [];
for i = 1:size(longscrewnodes,1)
    screwnodes = longscrewnodes(i,find(longscrewnodes(i,:)~=0));
    screwlength = norm(rncrystal(screwnodes(1),1:3)-rncrystal(screwnodes(end),1:3));
    if screwlength>L_cross_crit
        nodelist = [nodelist,screwnodes];
    end
end

if ~isempty(nodelist)
    L1=size(nodelist,2);
    [L2,L3]=size(connectivity);
    conlistall=zeros(L2,(L3-1)/2+1);
    conlistall(:,1)=connectivity(:,1);
    for i=1:L2
        connumb=conlistall(i,1);
        conlistall(i,2:connumb+1)=linspace(1,connumb,connumb);
    end
    
    conlist = conlistall(nodelist',:);
    
    for n=1:L1
        
        %                 if rncrystal(n,4) == 7 || rncrystal(n,4) == 67 || rncrystal(n,4) == 10 || rncrystal(n,4) == 6
        %                     continue
        %                 end
        
        if rncrystal(n,4) == 7 || rncrystal(n,4) == 67 || rncrystal(n,4) == 10
            continue
        end
        
        n0=nodelist(n);
        if rncrystal(n0,4) == 7
            continue
        end
        numNbrs=conlist(n,1);
        if numNbrs==2
            
            for i=1:numNbrs
                ii=conlist(n,i+1);
                linkid=connectivity(n0,2*ii);
                posinlink=connectivity(n0,2*ii+1);
                if posinlink == 1
                    n2=linkscrystal(linkid,3-posinlink);
                    segid2 = linkid;
                    posinlink2 = posinlink;
                    fsegn02=fsegcrystal(linkid,3*(posinlink-1)+[1:3]);
                    burgv2 = linkscrystal(linkid,3:5);
                    nplane2 = linkscrystal(linkid,6:8);
                else
                    n1=linkscrystal(linkid,3-posinlink);
                    segid1 = linkid;
                    posinlink1 = posinlink;
                    fsegn01=fsegcrystal(linkid,3*(posinlink-1)+[1:3]);
                    burgv1 = linkscrystal(linkid,3:5);
                    nplane1 = linkscrystal(linkid,6:8);
                end
            end
            
            reset_pos_n1 = 0;
            reset_pos_n2 = 0;
            reset_pos_n = 0;
            reset_plane_seg1 = 0;
            reset_plane_seg2 = 0;
            
            if ~exist('fsegn01','var') || ~exist('fsegn02','var')
                continue
            end
            
            fn0 = fsegn01+fsegn02;
            vec1 = rncrystal(n1,1:3)-rncrystal(n2,1:3);
            vec2 = rncrystal(n0,1:3)-rncrystal(n1,1:3);
            vec3 = rncrystal(n0,1:3)-rncrystal(n2,1:3);
            
            if norm(burgv1-burgv2)>eps
                disp('burgv is not conserved, see cross_slip_BCC')
                continue
            end
            burgv = burgv1;
            burg_mag = norm(burgv);
            burgv_norm = burgv./burg_mag;
            
            if ~(abs(abs(burgv(1))-abs(burgv(2))) < eps && abs(abs(burgv(1))-abs(burgv(3))) < eps &&...
                    abs(abs(burgv(2))-abs(burgv(3))) < eps)
                continue
            end
            
            test1 = dot(vec1,burgv_norm)^2;
            test2 = dot(vec2,burgv_norm)^2;
            test3 = dot(vec3,burgv_norm)^2;
            
            testmax1 = dot(vec1,vec1);
            testmax2 = dot(vec2,vec2);
            testmax3 = dot(vec3,vec3);
            
            seg1_is_screw = ((testmax2 - test2) < (testmax2 * s2thetacrit));
            seg2_is_screw = ((testmax3 - test3) < (testmax3 * s2thetacrit));
            bothseg_are_screw = (((testmax2 - test2) < (4.0 * testmax2 * s2thetacrit)) & ...
                ((testmax3 - test3) < (4.0 * testmax3 * s2thetacrit)) & ...
                ((testmax1 - test1) < (testmax1 * s2thetacrit)));
            
            if (seg1_is_screw || seg2_is_screw || bothseg_are_screw)
                Lseg1 = sqrt(testmax2);
                Lseg2 = sqrt(testmax3);
                fnodeThreshold = noiseFactor * MU * burg_mag * 0.5 * (Lseg1 + Lseg2);
                zipperThreshold1 = noiseFactor * MU * burg_mag *  Lseg1;
                zipperThreshold2 = noiseFactor * MU * burg_mag *  Lseg2;
                
                glideDir = eye(3)-burgv_norm.*burgv_norm';
                %             temp1 = glideDir.*nplane1;
                %             temp2 = glideDir.*nplane2;
                %             tempf = glideDir.*fn0;
                
                temp1 = glideDir*nplane1';
                temp2 = glideDir*nplane2';
                tempf = glideDir*fn0';
                
                plane1 = 1;
                plane2 = 1;
                fplane = 1;
                for j = 2:3
                    if abs(temp1(j)) < abs(temp1(plane1))
                        plane1 = j;
                    end
                    if abs(temp2(j)) < abs(temp2(plane2))
                        plane2 = j;
                    end
                    if abs(tempf(j)) > abs(tempf(fplane))
                        fplane = j;
                    end
                end
                
                newplane = cross(burgv_norm,glideDir(fplane,:));
                newplane = newplane ./ norm(newplane);
                
                if (bothseg_are_screw) & (plane1 == plane2) &...
                        (plane1 ~= fplane) &...
                        (abs(tempf(fplane)) > (weightFactor*abs(tempf(plane1))+fnodeThreshold))
                    
                    
                    
                    reset_plane_seg1 = 1;
                    reset_plane_seg2 = 1;
                    if ~(rncrystal(n1,4) == 7 || rncrystal(n1,4) == 10)
                        new_pos_n1 = rncrystal(n2,1:3) + dot(vec1,burgv_norm).*burgv_norm;
                        reset_pos_n1 = 1;
                        new_pos_n = rncrystal(n2,1:3) + dot(vec3,burgv_norm).*burgv_norm;
                        reset_pos_n = 1;
                    elseif ~(rncrystal(n2,4) == 7 || rncrystal(n2,4) == 10)
                        new_pos_n2 = rncrystal(n1,1:3) - dot(vec1,burgv_norm).*burgv_norm;
                        reset_pos_n2 = 1;
                        new_pos_n = rncrystal(n1,1:3) + dot(vec2,burgv_norm).*burgv_norm;
                        reset_pos_n = 1;
                    end
                    fdotglide = dot(fn0,glideDir(fplane,:));
                    temp = areamin / abs(dot(vec1,burgv_norm)) * 2.0 * (1.0 + eps) * sign(fdotglide);
                    new_pos_n = new_pos_n + 1.*temp.*glideDir(fplane,:);
                    rntemp = rn;
                    rntemp(n0,1:3) = new_pos_n*rotationBCC;
                    if reset_pos_n1 == 1
                        rntemp(n1,1:3) = new_pos_n1*rotationBCC;
                    end
                    if reset_pos_n2 == 1
                        rntemp(n2,1:3) = new_pos_n2*rotationBCC;
                    end
                    linkstemp = links;
                    linkstemp(segid1,6:8) = newplane*rotationBCC;
                    linkstemp(segid2,6:8) = newplane*rotationBCC;
                    fsegtemp = segforcevec(MU,NU,a,Ec,rntemp,linkstemp,0,vertices,...
    uhat,nc,xnodes,D,mx,mz,w,h,d);
                    fsegnew(:,1:3) = fsegtemp(:,1:3)*rotationBCC';
                    fsegnew(:,4:6) = fsegtemp(:,4:6)*rotationBCC';
                    fn0new = fsegnew(segid1,3*(posinlink1-1)+[1:3])+fsegnew(segid2,3*(posinlink2-1)+[1:3]);
                    fdotglidenew = dot(fn0new,glideDir(fplane,:));
                    if sign(fdotglidenew) * sign(fdotglide) < 0.0
                        reset_pos_n1 = 0;
                        reset_pos_n2 = 0;
                        reset_pos_n = 0;
                        reset_plane_seg1 = 0;
                        reset_plane_seg2 = 0;
                    end
                elseif (seg1_is_screw) & (plane1 ~= plane2) &...
                        (plane2 == fplane) & ...
                        (abs(tempf(fplane)) > (weightFactor*abs(tempf(plane1))+fnodeThreshold))
                    %             elseif (seg1_is_screw) & (plane1 ~= plane2) &...
                    %                     (plane2 == fplane)
                    fdotplane1 = abs(glideDir(plane1,:)*fsegn01');
                    fdotplanef = abs(glideDir(fplane,:)*fsegn01');
                    if (fdotplanef > zipperThreshold1 + fdotplane1)
                        reset_pos_n1 = 1;
                        reset_plane_seg1 = 1;
                        new_pos_n1 = rncrystal(n0,1:3) - dot(vec2,burgv_norm).*burgv_norm;
                    end
                elseif (seg2_is_screw) && (plane1 ~= plane2) &...
                        (plane1 == fplane) & ...
                        (abs(tempf(fplane)) > (weightFactor*abs(tempf(plane2))+fnodeThreshold))
                    %             elseif (seg2_is_screw) & (plane1 ~= plane2) &...
                    %                     (plane1 == fplane)
                    fdotplane2 = abs(glideDir(plane2,:)*fsegn02');
                    fdotplanef = abs(glideDir(fplane,:)*fsegn02');
                    if (fdotplanef > zipperThreshold2 + fdotplane2)
                        reset_pos_n2 = 1;
                        reset_plane_seg2 = 1;
                        new_pos_n2 = rncrystal(n0,1:3) - dot(vec3,burgv_norm).*burgv_norm;
                    end
                end
            end
            %                     if norm(rn(n,1:3))>-0.1
            if reset_pos_n1 == 1
                rn(n1,1:3) = new_pos_n1*rotationBCC;
            end
            if reset_pos_n == 1
                rn(n0,1:3) = new_pos_n*rotationBCC;
            end
            if reset_pos_n2 == 1
                rn(n2,1:3) = new_pos_n2*rotationBCC;
            end
            
            %                         if abs(newplane(3))>-0.1
            if reset_plane_seg1 == 1
                links(segid1,6:8) = newplane*rotationBCC;
                disp('seg1 cross slips on new plane:'); newplane
%                 do_reset = 0;
%                 last_step = curstep;
            end
            if reset_plane_seg2 == 1
                links(segid2,6:8) = newplane*rotationBCC;
                disp('seg2 cross slips on new plane:'); newplane
%                 do_reset = 0;
%                 last_step = curstep;
            end
            %                         end
            
            if (reset_plane_seg1 == 1 || reset_plane_seg2 == 1) & ~(reset_plane_seg1 == 1 & reset_plane_seg2 == 1)
                disp('zipper')
            end
            %                     end
        end
    end
end


% if reset_plane_seg1 == 1 || reset_plane_seg2 == 1
%     do_reset = 0;
%     last_step = curstep;
% end
% 
% if curstep-last_step>0
%     do_reset = 1;
% end

end
