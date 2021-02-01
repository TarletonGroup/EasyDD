%% Find surface nodes.
amag = 3.18e-4; 
surfnIndices = find(rn(:,4)==6);
surflIndices = [];
for i = 1:size(surfnIndices,1)
    tmp = find(links(:,1)==surfnIndices(i));
    surflIndices = [surflIndices; tmp];
end
rn2 = rn;
rn2(:,1:3) = amag*rn2(:,1:3); 
surfNodes = rn2(surfnIndices,:);
surfLinks = links(surflIndices,:);

condx = (surfNodes(:,1) > 9.9);
surfNodes2 = surfNodes(condx,:);

nodelabel = find(rn2(:,1)==surfNodes2(1));

test = rn2(find(rn2(:,1)==surfNodes2(1)),:);

linksnode = [links(links(:,1) == nodelabel,:); links(links(:,2) == nodelabel,:)];
linksnode(:,1:2)

rn2(88,:)
rn2(1650,:)

