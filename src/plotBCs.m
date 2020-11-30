function plotBCs(saveBC, curstep)
    %===============================================================%
    % Daniel Hortelano Roig (11/11/2020)
    % daniel.hortelanoroig@materials.ox.ac.uk 
    
    % Plots relationship between BCs as time-steps progress.
    %===============================================================%
    
    %% Extraction
    
    Xsim = saveBC.Xsim; Ysim = saveBC.Ysim;
    Xlabel = saveBC.Xlabel; Ylabel = saveBC.Ylabel;
    figname = saveBC.figname;
    
    %% Plotting
    
    % If an axis is 0-dimension, don't plot:
    if (min(Xsim) == max(Xsim)) || (min(Ysim) == max(Ysim))
        return
    end
    
    % Setup figure:
    figure('Name',figname);
    grid on, hold on
    
    % Plotting:
    plot([min(Xsim) max(Xsim)], zeros(2),'--','LineWidth',1,'Color','black') % x-axis
    plot(zeros(2), [min(Ysim) max(Ysim)],'--','LineWidth',1,'Color','black') % y-axis
    plot(Xsim(1:curstep), Ysim(1:curstep), ... % Main plot
        '-','LineWidth',3,'Color','magenta')
    
    % Axes and labels:
    ax = gca; % Define axes object
    ax.XLim = [min(Xsim) max(Xsim)]; % Define limits on axes
    ax.YLim = [min(Ysim) max(Ysim)];
    xlabel(Xlabel,'FontSize',14); % Label axes
    ylabel(Ylabel,'FontSize',14);
    title(figname,'FontSize', 18); % Plot title
    legend('x-axis','y-axis','BCs'); % Legend
    
end
