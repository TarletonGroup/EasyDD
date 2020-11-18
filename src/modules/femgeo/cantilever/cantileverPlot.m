function cantileverPlot(BCstore, curstep, simscale)
    %===============================================================%
    % Daniel Hortelano Roig (11/11/2020)
    % daniel.hortelanoroig@materials.ox.ac.uk 
    
    % Plots relationship between displacement and force BCs as
    % the time-steps progress. The scale structure stores the BCs
    % and their scaling for the plot.
    %===============================================================%
    
    %% Extraction
    
    Usim = BCstore.Usim; Fsim = BCstore.Fsim;
    
    %% Plotting
    
    % Setup figure:
    cantileverPlot_figname = "Bending cantilever BCs";
    figure('Name',cantileverPlot_figname);
    grid on, hold on
    
    % Plotting:
    plot(Usim(1:curstep), -Fsim(1:curstep), ... % Main plot
        '-','LineWidth',3,'Color','magenta')
    plot([min(Usim) max(Usim)], zeros(2),'--','LineWidth',1,'Color','black') % x-axis
    plot(zeros(2), [min(Fsim) max(Fsim)],'--','LineWidth',1,'Color','black') % y-axis
    
    % Axes and labels:
    ax = gca; % Define axes object
    ax.XLim = [min(Usim) max(Usim)]; % Define limits on axes
    ax.YLim = [min(Fsim) max(Fsim)];
    xlabel('Force','FontSize',14); % Label axes
    ylabel('Displacement','FontSize',14);
    title(cantileverPlot_figname,'FontSize',20); % Title plot
    legend('x-axis','y-axis','8b'); % Legend
    
end
