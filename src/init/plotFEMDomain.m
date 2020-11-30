function plotFEMDomain(S, FEM)
    %===============================================================%
    % Daniel Hortelano Roig (29/11/2020)
    % daniel.hortelanoroig@materials.ox.ac.uk 
    
    % Constructs plot of the geometry before the simulation begins.
    %===============================================================%
    
    %% Extraction
    
    % Surface sets:
    Scat = S.cat;
    
    % FEM:
    xnodes = FEM.xnodes;
    
    %% Plot
    
    % Setup figure:
    figname = "Initial domain surface nodes";
    figure('Name',figname);
    grid on, hold on
    
    % Plotting:
    plot3([min(xnodes(Scat(:,1),1)) max(xnodes(Scat(:,1),1))], zeros(2),zeros(2), ... % x-axis
        '--','LineWidth',1,'Color','black')
    plot3(zeros(2), [min(xnodes(Scat(:,1),2)) max(xnodes(Scat(:,1),2))], zeros(2), ... % y-axis
        '--','LineWidth',1,'Color','black')
    plot3(zeros(2), zeros(2), [min(xnodes(Scat(:,1),3)) max(xnodes(Scat(:,1),3))], ... % z-axis
        '--','LineWidth',1,'Color','black')
    plot3(xnodes(Scat(:, 1), 1), xnodes(Scat(:, 1), 2), xnodes(Scat(:, 1), 3), ... % Main plot
        '*', 'MarkerSize',1,'Color','red')
    
    % Axes and labels:
    xlabel('x','FontSize',14); % Label axes
    ylabel('y','FontSize',14);
    zlabel('z','FontSize',14);
    title(figname,'FontSize',18); % Plot title
    
    axis('equal')
end
