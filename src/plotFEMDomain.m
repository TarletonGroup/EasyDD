function fig = plotFEMDomain(Stop, Sbot, Sright, Sleft, Sfront, Sback, Smixed, xnodes)
    fig = figure(1); clf; hold on; view(3)
    xlabel('x'); ylabel('y'); zlabel('z')
    plot3(xnodes(Stop(:, 1), 1), xnodes(Stop(:, 1), 2), xnodes(Stop(:, 1), 3), 'r*')
    plot3(xnodes(Sbot(:, 1), 1), xnodes(Sbot(:, 1), 2), xnodes(Sbot(:, 1), 3), 'r*')
    plot3(xnodes(Sright(:, 1), 1), xnodes(Sright(:, 1), 2), xnodes(Sright(:, 1), 3), 'b.')
    plot3(xnodes(Sleft(:, 1), 1), xnodes(Sleft(:, 1), 2), xnodes(Sleft(:, 1), 3), 'b.')
    plot3(xnodes(Sfront(:, 1), 1), xnodes(Sfront(:, 1), 2), xnodes(Sfront(:, 1), 3), 'k*')
    plot3(xnodes(Sback(:, 1), 1), xnodes(Sback(:, 1), 2), xnodes(Sback(:, 1), 3), 'k*')
    plot3(xnodes(Smixed(:, 1), 1), xnodes(Smixed(:, 1), 2), xnodes(Smixed(:, 1), 3), 'g*')
    axis('equal')
    hold off
end
