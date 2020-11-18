function fig = plotFEMDomain(S, xnodes)
    
    fig = figure(1); clf; hold on; view(3)
    xlabel('x'); ylabel('y'); zlabel('z')
    plot3(xnodes(S.top(:, 1), 1), xnodes(S.top(:, 1), 2), xnodes(S.top(:, 1), 3), 'r*')
    plot3(xnodes(S.bot(:, 1), 1), xnodes(S.bot(:, 1), 2), xnodes(S.bot(:, 1), 3), 'r*')
    plot3(xnodes(S.right(:, 1), 1), xnodes(S.right(:, 1), 2), xnodes(S.right(:, 1), 3), 'b.')
    plot3(xnodes(S.left(:, 1), 1), xnodes(S.left(:, 1), 2), xnodes(S.left(:, 1), 3), 'b.')
    plot3(xnodes(S.front(:, 1), 1), xnodes(S.front(:, 1), 2), xnodes(S.front(:, 1), 3), 'k*')
    plot3(xnodes(S.back(:, 1), 1), xnodes(S.back(:, 1), 2), xnodes(S.back(:, 1), 3), 'k*')
    plot3(xnodes(S.topright(:, 1), 1), xnodes(S.topright(:, 1), 2), xnodes(S.topright(:, 1), 3), 'g*')
    axis('equal')
    hold off
end
