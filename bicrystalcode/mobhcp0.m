function [vn,fn] = mobhcp0(fseg,rn,links,connectivity,nodelist,conlist)
% mobility law function (model: HCP)


global HCP_A_Basal_EdgeDrag HCP_A_Basal_ScrewDrag HCP_A_Prismatic_EdgeDrag HCP_A_Prismatic_ScrewDrag HCP_A_1stPyramidal_EdgeDrag
global HCP_A_1stPyramidal_ScrewDrag HCP_A_2ndPyramidal_EdgeDrag HCP_A_2ndPyramidal_ScrewDrag HCP_CpA_Prismatic_EdgeDrag
global HCP_CpA_Prismatic_ScrewDrag HCP_CpA_1stPyramidal_EdgeDrag HCP_CpA_1stPyramidal_ScrewDrag HCP_CpA_2ndPyramidal_EdgeDrag
global HCP_CpA_2ndPyramidal_ScrewDrag HCP_C_Prismatic_EdgeDrag HCP_C_Prismatic_ScrewDrag HCP_Sessile_EdgeDrag
global HCP_Sessile_ScrewDrag HCP_LineDrag

global cOVERa
global MU NU maxSeg
global burgsref planesref edgesref bpiref ppbref planestyperef



% Parameter for line mobility: fast mobility, small drag
dragLine = HCP_LineDrag;


% Parameter for climb mobility: slow mobility, large drag
dragClimb = HCP_Sessile_ScrewDrag;

% Set up some arrays with the various drag coefficients from the
% control parameters in order to get the proper plane type
% later.

AdragScrew(1) = HCP_A_Basal_ScrewDrag;
AdragScrew(2) = HCP_A_Prismatic_ScrewDrag;
AdragScrew(3) = HCP_A_1stPyramidal_ScrewDrag;
AdragScrew(4) = HCP_A_2ndPyramidal_ScrewDrag;
AdragScrew(5) = HCP_Sessile_ScrewDrag;

CdragScrew(1) = HCP_C_Prismatic_ScrewDrag;
CdragScrew(2) = HCP_Sessile_ScrewDrag;

CpAdragScrew(1) = HCP_CpA_Prismatic_ScrewDrag;
CpAdragScrew(2) = HCP_CpA_1stPyramidal_ScrewDrag;
CpAdragScrew(3) = HCP_CpA_2ndPyramidal_ScrewDrag;
CpAdragScrew(4) = HCP_Sessile_ScrewDrag;

AdragEdge(1) = HCP_A_Basal_EdgeDrag;
AdragEdge(2) = HCP_A_Prismatic_EdgeDrag;
AdragEdge(3) = HCP_A_1stPyramidal_EdgeDrag;
AdragEdge(4) = HCP_A_2ndPyramidal_EdgeDrag;
AdragEdge(5) = HCP_Sessile_EdgeDrag;

CdragEdge(1) = HCP_C_Prismatic_EdgeDrag;
CdragEdge(2) = HCP_Sessile_EdgeDrag;

CpAdragEdge(1) = HCP_CpA_Prismatic_EdgeDrag;
CpAdragEdge(2) = HCP_CpA_1stPyramidal_EdgeDrag;
CpAdragEdge(3) = HCP_CpA_2ndPyramidal_EdgeDrag;
CpAdragEdge(4) = HCP_Sessile_EdgeDrag;
        

% numerical tolerance
eps=1e-12;

% length of the nodelist for which the velocity will be calculated
L1=size(nodelist,1);

% if no nodelist is given then the nodelist becomes the whole node population
% this portion of the code sets up that nodelist along with the connlist
% that contains all of the nodal connections

if L1==0
    L1=size(rn,1);
    nodelist=linspace(1,L1,L1)';
    [L2,L3]=size(connectivity);
    conlist=zeros(L2,(L3-1)/2+1);
    conlist(:,1)=connectivity(:,1);
    for i=1:L2
        connumb=conlist(i,1);
        conlist(i,2:connumb+1)=linspace(1,connumb,connumb);
    end
end 

% Loop over all the nodes
for n=1:L1
    % n0 is the nodeid of the nth node in nodelist
    n0=nodelist(n);                 
    % numNbrs is the number of connections for node n0 in conlist
    numNbrs=conlist(n,1);           
    % initialize the total force and the total drag matrix
    fn(n,:)=zeros(1,3);             

    linedir = zeros(numNbrs,3);
    burg    = zeros(numNbrs,3);
    junctionExists = zeros(numNbrs,1);
    edgeExists     = zeros(numNbrs,1);
    screwExists    = zeros(numNbrs,1);

    dragGlideScrew = zeros(numNbrs,1);     
    dragGlideEdge1 = zeros(numNbrs,1);     
    dragGlideEdge2 = zeros(numNbrs,1);     

    edgeDir1 = zeros(numNbrs,3);
    edgeDir2 = zeros(numNbrs,3);

    segScrewLength = zeros(numNbrs,1); 
    segEdge1Length = zeros(numNbrs,1); 
    segEdge2Length = zeros(numNbrs,1); 

    ndir = zeros(numNbrs,3);
    segLen = zeros(numNbrs,1);

    for i=1:numNbrs
        ii=conlist(n,i+1);                                                                      
        % connectionid for this connection
        linkid=connectivity(n0,2*ii);
        posinlink=connectivity(n0,2*ii+1);
        n1=links(linkid,3-posinlink);

        % calculate the length of the link and its tangent line direction
        rt=rn(n1,1:3)-rn(n0,1:3);                                                               
        segLen(i)=norm(rt);
        if segLen(i)>0.0
            fsegn0=fseg(linkid,3*(posinlink-1)+[1:3]);
            fn(n,:)=fn(n,:)+fsegn0; 

            burg(i,:)    = links(connectivity(n0,2*ii),3:5); 
            plane   = links(connectivity(n0,2*ii),6:8); 
            linedir(i,:) = rt./segLen(i);


            % find which Burgers vector is the current one in the reference
            % list.
            bIndex = -1;
            for k=1:length(burgsref)
                sbsign = dot(burg(i,:),burgsref(k,:));                
                sb = sign(sbsign);

                if (abs(norm(burg(i,:)) - norm(burgsref(k,:))) < 1e-4 && ...
                    abs(burg(i,1)-sb*burgsref(k,1))<1e-4 && ...
                    abs(burg(i,2)-sb*burgsref(k,2))<1e-4 && ...
                    abs(burg(i,3)-sb*burgsref(k,3))<1e-4)
                     bIndex = k;
                     break;
                end
            end

            % find the list of glide planes
            if (bIndex > 0) 
                % All possible planes for that Burgers vector
                glideplanes = planesref(bpiref(bIndex):bpiref(bIndex)+ppbref(bIndex)-1,:);

                % Careful : pIndex is a value between -1 and 4;
                pIndex = -1;
                for k=1:size(glideplanes,1)
                    spsign = dot(plane,glideplanes(k,:));                
                    sp = sign(spsign);

                    if (abs(norm(plane) - norm(glideplanes(k,:))) < 1e-4 && ...
                        abs(plane(1)-sp*glideplanes(k,1))<1e-4 && ...
                        abs(plane(2)-sp*glideplanes(k,2))<1e-4 && ...
                        abs(plane(3)-sp*glideplanes(k,3))<1e-4)
                         pIndex = k;
                         break;
                    end
                end
                
                if (pIndex > 0)
                    ndir(i,:) = glideplanes(pIndex,:);
                    planeType = planestyperef(bpiref(bIndex)-1+pIndex);
                else
                    planeType = 5;
                end
            else
                pIndex = -1;
                planeType = 5;
            end

            if (bIndex <= 0) 
                % Junction
                junctionExists(i) = 1;
                edgeExists(i)     = 0;
                screwExists(i)    = 0;
            elseif (planeType == 5)
                % Sessile plane : treated as a junction
                junctionExists(i) = 1;
                edgeExists(i)     = 0;
                screwExists(i)    = 0;
            else
               junctionExists(i) = 0;
               
               % Assign screw glide
               if (bIndex <= 3)
                   dragGlideScrew(i) = AdragScrew(planeType);
               elseif (bIndex == 10)
                   if (planeType == 2)
                    dragGlideScrew(i) = CdragScrew(1);
                   else  
                    dragGlideScrew(i) = CdragScrew(2);
                   end
               else
                  dragGlideScrew(i) = CpAdragScrew(planeType-1); 
               end

               % Assign edge glide
               edgeIndex1 = bpiref(bIndex)-1 + pIndex; 

               [segScrewLength(i), segEdge1Length(i), segEdge2Length(i), edgeDir1(i,:), edgeDir2(i,:), edgeIndex2] = decomposearm(rt,bIndex,burgsref,edgesref,bpiref,ppbref,pIndex);

               planeType1 = planestyperef(edgeIndex1);
               planeType2 = planestyperef(edgeIndex2);

               if (bIndex <= 3)
                   dragGlideEdge1(i) = AdragEdge(planeType1);
                   dragGlideEdge2(i) = AdragEdge(planeType2);
                elseif (bIndex == 10)
                   if (planeType == 2)
                    dragGlideEdge1(i) = CdragEdge(1);
                    dragGlideEdge2(i) = CdragEdge(1);
                   else  
                     dragGlideEdge1(i) = CdragEdge(2);
                     dragGlideEdge2(i) = CdragEdge(2);
                   end
               else
                  dragGlideEdge1(i) = CpAdragEdge(planeType1-1); 
                  dragGlideEdge2(i) = CpAdragEdge(planeType2-1); 
               end

                edgeLength = sqrt((segEdge1Length(i) * segEdge1Length(i)) + (segEdge2Length(i) * segEdge2Length(i)));
                
                edgeExists(i)  = ((edgeLength / segLen(i)) > 1.0e-04);
                screwExists(i) = ((segScrewLength(i)/segLen(i))>1.0e-04);
            end % bIndex 
        end      
    end % i = 1:numNbrs


  % START NEWTON-RAPHSON
    count = 0;
    maxiter = 5;

    vtest = zeros(1,3);
    ferror=fn(n,:);
    dferrordv=zeros(3,3);
    forceerror=(ferror*ferror');
    lmax = maxSeg;
    ferrortol = max((fn(n,:)*fn(n,:)'*1e-16),(MU*MU*lmax*lmax*1e-24));

    while ( (count < maxiter) & (forceerror > ferrortol) )
        count = count + 1;

        ferror=fn(n,:);
        dferrordv=zeros(3,3);


        for i=1:numNbrs   
            if (segLen(i) > 0)                 
                % junction contribution
                if (junctionExists(i))
                    [fjunc, dfjuncdv] = JunctionDrag(vtest, linedir(i,:), dragLine, dragClimb);
                    ferror=ferror - (0.5*segLen(i)).*fjunc;
                    dferrordv=dferrordv - (0.5*segLen(i)).*dfjuncdv;
                else
                    % screw contribution
                    if (screwExists(i))                            
                        [fScrew, dfScrewdv] = ScrewDrag(vtest, burg(i,:), ndir(i,:), dragClimb, dragLine, dragGlideScrew(i));
                        ferror = ferror - 0.5 * segScrewLength(i) * fScrew;
                        dferrordv = dferrordv - 0.5 * segScrewLength(i)*dfScrewdv;
                    end

                    % edge contribution 
                    if (edgeExists(i))                   
                        % for first edge 
                        [fEdge, dfEdgedv] = EdgeDrag(vtest, burg(i,:), edgeDir1(i,:),dragClimb, dragLine, dragGlideEdge1(i,:));
                        ferror = ferror - 0.5 * segEdge1Length(i) * fEdge;
                        dferrordv = dferrordv - 0.5 *segEdge1Length(i) *dfEdgedv;  

                        % for second edge 
                        [fEdge, dfEdgedv] = EdgeDrag(vtest, burg(i,:), edgeDir2(i,:),dragClimb, dragLine, dragGlideEdge2(i,:));
                        ferror = ferror - 0.5 * segEdge2Length(i) * fEdge;
                        dferrordv = dferrordv - 0.5 *segEdge2Length(i) *dfEdgedv;
                    end
                end
            end
        end % loop over segments

        if (abs(det(dferrordv)) < eps) 
            count
            burg(i,:)
            segLen(i)
            segScrewLength(i)
            segEdge1Length(i)
            segEdge2Length(i) 
            ferror
            screwExists(i)
            edgeExists(i)
            dferrordv
            det(dferrordv)
            pause
        end
        vtest=vtest-(dferrordv\ferror')';
        forceerror =  (ferror*ferror');

        %fprintf('count=%d forceerror=%15.10e ferrortol=%15.10e det(dferrordv)=%15.10e\n',count,forceerror,ferrortol,det(dferrordv));
    end

    vn(n,:)=vtest;

    % if Newton-Raphson has not converged for the node n, report a
    % failed to find a velocity

   if ( (count >= maxiter) || (forceerror > ferrortol) )
    fprintf('Newton-Raphson has failed\n');
    fprintf('count=%d det(dferrordv)=%f  forceerror=%f ferrortol=%f\n',count,det(dferrordv),forceerror,ferrortol);
    pause
   end

   end  % end loop over all ns
end


function [fEdgeDrag, dfEdgeDragdv] = EdgeDrag(vel, burg, edgeDir, dragClimb, dragLine, dragGlide)
    glideDir = burg/norm(burg);
    lineDir  = edgeDir/norm(edgeDir);
    climbDir = cross(lineDir, glideDir);
    climbDir = climbDir/norm(climbDir);
    
    glideDirMat = glideDir' * glideDir;
    climbDirMat = climbDir' * climbDir;
    lineDirMat  = lineDir'  * lineDir;
    
    dfEdgeDragdv = dragClimb * climbDirMat + dragLine  * lineDirMat  + dragGlide * glideDirMat;
    fEdgeDrag = (dfEdgeDragdv * vel')';
end



function [fScrewDrag, dfScrewDragdv] = ScrewDrag(vel, burg, climbDir, dragClimb, dragLine, dragGlide)
    lineDir  = burg/norm(burg);
    glideDir = cross(climbDir, lineDir);
    glideDir = glideDir/norm(glideDir);

    glideDirMat = glideDir' * glideDir;
    climbDirMat = climbDir' * climbDir;
    lineDirMat  = lineDir'  * lineDir;
        
    dfScrewDragdv = dragClimb * climbDirMat + dragLine  * lineDirMat  + dragGlide * glideDirMat; 
    fScrewDrag = (dfScrewDragdv * vel')';
end


function [fJunctionDrag, dfJunctionDragdv] = JunctionDrag(vel, lineDir, BJunctionLine, BJunctionGlide)
    glideDirMat = ones(3,3) - lineDir' * lineDir; 
    lineDirMat  = lineDir'  * lineDir;

    dfJunctionDragdv = BJunctionGlide * glideDirMat + BJunctionLine  * lineDirMat;
    fJunctionDrag    = (dfJunctionDragdv * vel')';
end

function [ls,le1,le2,te1,te2,e2] = decomposearm(seg,bIndex,burgref,edgeref,bpi,ppb,alpha)

% decomposition of seg 
bn = burgref(bIndex,:);
bn = bn/norm(bn);

% screw 
ls = dot(bn,seg);

% edges
segleft = seg - ls*bn;

% get indices (alpha, beta) of the two highest values
edge = edgeref(bpi(bIndex):bpi(bIndex)+ppb(bIndex)-1,1:3); % 3 or 4 edges

ltemp = zeros(ppb(bIndex),1);
for i=1:ppb(bIndex)
    ltemp(i) = abs(dot(edge(i,:),segleft'));
end

%[lmax,alpha] = max(ltemp);
beta = Get2ndMaxIndex(ltemp,alpha);

te1 = edge(alpha,:);
te2 = edge(beta,:);

[le1,le2] = DecompVec(segleft,te1,te2);

% Positive lengths
le1 = abs(le1);
le2 = abs(le2);

% indices of the planes associated with the hightest edges.
% e1= bpi(bIndex)+alpha-1;
e2= bpi(bIndex)+beta-1;
end


function [indexOf2ndMax]=Get2ndMaxIndex(ltemp,index)
    indexOf2ndMax = (index == 1) + 1;
        
    for m=1:length(ltemp)                
        if ( (ltemp(m) > ltemp(indexOf2ndMax)) && (m ~= index) )
            indexOf2ndMax = m;
        end
    end
end


function [outV1,outV2] = DecompVec(inVec, vec1, vec2)
    cc   = dot(vec1,  vec2);
    dec1 = dot(inVec, vec1);
    dec2 = dot(inVec, vec2);
    
    outV1 = (dec1 - (dec2 * cc)) / (1.0 - (cc * cc));
    outV2 = (dec2 - (dec1 * cc)) / (1.0 - (cc * cc));
end































