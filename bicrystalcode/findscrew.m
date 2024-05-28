function [longscrewsegs,longscrewnodes] = findscrew(rn,links,connectivity,s2thetacrit)

longscrewsegs = [];
longscrewnodes = [];

L1=size(links,1);

k=0;
screwseg =[];
screw_linked = [];
for i=1:L1
    if connectivity(links(i,1),1)>2 || connectivity(links(i,2),1)>2
        continue
    end%HY20191016: jump over junctions
    ldir = rn(links(i,1),1:3)-rn(links(i,2),1:3);
    ldir = ldir/norm(ldir);
    bvec = links(i,3:5);
    bvec = bvec/norm(bvec);
    if 1-dot(ldir,bvec)^2<s2thetacrit
        k = k+1;
        screwseg(k) = i;
    end
end

if ~isempty(screwseg)
    L2 = size(screwseg,2);
    n_screw_linked = 0;
    for i=1:L2
        n1 = links(screwseg(i),1);
        n2 = links(screwseg(i),2);
        for j = 1:L2
            n1j = links(screwseg(j),1);
            n2j = links(screwseg(j),2);
            if (n1==n1j || n1==n2j || n2==n1j || n2==n2j) && (i~=j)
                n_screw_linked = n_screw_linked+1;
                screw_linked(n_screw_linked) = screwseg(i);
                break
            end
        end
    end
    
    if ~isempty(screw_linked)
        L3 = size(screw_linked,2);
        n_endnodes = 0;
        for i=1:L3
            allnodes(i*2-1) = links(screw_linked(i),1);
            allnodes(i*2) = links(screw_linked(i),2);
        end
        
        for i=1:2*L3
            ntemp = find(allnodes==allnodes(i));
            if size(ntemp,2)==1
                n_endnodes = n_endnodes+1;
                endnodes(n_endnodes,1) = allnodes(i);
                endnodes(n_endnodes,2) = 0; %HY20191016: a flag
            end
        end
        
        L4 = size(endnodes,1);
        n_longscrew = 0;
        for i=1:L4
            if endnodes(i,2) == 0
                endnodes(i,2) = 1;
                n_longscrew = n_longscrew+1;
                targetnode = endnodes(i,1);
                longscrewnodes(n_longscrew,1) = targetnode;
                targetsegid = 0;
                nn = 0;
                next = 1;
                while next == 1
                    found = 0;
                    for j = 1:L3
                        if targetnode == links(screw_linked(j),1) && targetsegid~=screw_linked(j)
                            targetnode = links(screw_linked(j),2);
                            targetsegid = screw_linked(j);
                            found = 1;
                            break
                        elseif targetnode == links(screw_linked(j),2) && targetsegid~=screw_linked(j)
                            targetnode = links(screw_linked(j),1);
                            targetsegid = screw_linked(j);
                            found = 1;
                            break
                        end
                    end
                    if found == 1
                        nn = nn+1;
                        longscrewsegs(n_longscrew,nn) = targetsegid;
                        longscrewnodes(n_longscrew,nn+1) = targetnode;
                    else
                        next = 0;
                    end
                end
                endnodes(find(endnodes(:,1)==targetnode),2) = 1;
            end
        end
    end
end

% common = intersect(Sright(:,1),Scontact(:,1));
% Sright(ismember(Sright(:,1), common),:) = [];

end