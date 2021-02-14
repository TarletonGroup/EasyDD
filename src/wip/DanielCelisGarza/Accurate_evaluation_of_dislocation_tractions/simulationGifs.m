clear all

nameAIni = '../../../../output/initial_analytic_11-Jan-2021_bb_0';
nameA = '../../../../output/analytic_11-Jan-2021_bb_';
nameNIni = '../../../../output/initial_numeric_11-Jan-2021_bb_0';
nameN = '../../../../output/numeric_11-Jan-2021_bb_';
addpath '../../../'

%%
axis tight manual % this ensures that getframe() returns a consistent size
load(nameAIni);
counter = 0;
outputname = './paper/images/analytic.gif';
while counter < 66
    filename = strcat(nameA, sprintf('%d', counter));
    load(filename)
    view([0 0])
%     view([-15 15])
    plotnodes(rn,links,plim,vertices);
    % Capture the plot as an image
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if counter == 0
        imwrite(imind,cm,outputname,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,outputname,'gif','WriteMode','append');
    end
    counter = counter + 1;
end

%%
axis tight manual
load(nameAIni)
counter = 0;
outputname = './paper/images/analytic_';
filename = strcat(nameA, sprintf('%d', counter));
outputname = strcat(outputname, sprintf('%d', counter));
load(filename)
plotnodes(rn,links,plim,vertices);
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(gcf, sprintf('%s.pdf', outputname), '-dpdf', '-r300')

counter = 30;
outputname = './paper/images/analytic_';
filename = strcat(nameA, sprintf('%d', counter));
outputname = strcat(outputname, sprintf('%d', counter));
load(filename)
plotnodes(rn,links,plim,vertices);
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(gcf, sprintf('%s.pdf', outputname), '-dpdf', '-r300')

counter = 50;
outputname = './paper/images/analytic_';
filename = strcat(nameA, sprintf('%d', counter));
outputname = strcat(outputname, sprintf('%d', counter));
load(filename)
plotnodes(rn,links,plim,vertices);
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(gcf, sprintf('%s.pdf', outputname), '-dpdf', '-r300')

counter = 66;
outputname = './paper/images/analytic_';
filename = strcat(nameA, sprintf('%d', counter));
outputname = strcat(outputname, sprintf('%d', counter));
load(filename)
plotnodes(rn,links,plim,vertices);
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(gcf, sprintf('%s.pdf', outputname), '-dpdf', '-r300')

%%
axis tight manual % this ensures that getframe() returns a consistent size
load(nameNIni);
counter = 0;
outputname = './paper/images/numeric.gif';
while counter < 156
    filename = strcat(nameN, sprintf('%d', counter));
    load(filename)
    view([0 0])
%     view([-15 15])
    plotnodes(rn,links,plim,vertices);
    % Capture the plot as an image
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if counter == 0
        imwrite(imind,cm,outputname,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,outputname,'gif','WriteMode','append');
    end
    counter = counter + 2;
end

%%
axis tight manual
load(nameNIni)
counter = 0;
outputname = './paper/images/numeric_';
filename = strcat(nameN, sprintf('%d', counter));
outputname = strcat(outputname, sprintf('%d', counter));
load(filename)
plotnodes(rn,links,plim,vertices);
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(gcf, sprintf('%s.pdf', outputname), '-dpdf', '-r300')

counter = 78;
outputname = './paper/images/numeric_';
filename = strcat(nameN, sprintf('%d', counter));
outputname = strcat(outputname, sprintf('%d', counter));
load(filename)
plotnodes(rn,links,plim,vertices);
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(gcf, sprintf('%s.pdf', outputname), '-dpdf', '-r300')

counter = 118;
outputname = './paper/images/numeric_';
filename = strcat(nameN, sprintf('%d', counter));
outputname = strcat(outputname, sprintf('%d', counter));
load(filename)
plotnodes(rn,links,plim,vertices);
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(gcf, sprintf('%s.pdf', outputname), '-dpdf', '-r300')

counter = 156;
outputname = './paper/images/numeric_';
filename = strcat(nameN, sprintf('%d', counter));
outputname = strcat(outputname, sprintf('%d', counter));
load(filename)
plotnodes(rn,links,plim,vertices);
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])
print(gcf, sprintf('%s.pdf', outputname), '-dpdf', '-r300')
