function set_plot( current_fig, current_axes )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% To avoid figures with random size (matlab bug)
pause(0.1);    
drawnow; 

colors = [
    1.00  0.00  0.00
    0.00  0.50  0.00
    0.00  0.00  1.00
    0.93  0.69  0.13
    0.50  0.00  0.50
    ];

lineStyle = {'-','--',':','-.'};
assignin('caller','lineStyle',lineStyle);

set(current_fig, ...
    'Position',[100 100 800 600], ...
    'Color','white', ...
    'DefaultTextInterpreter', 'latex')

set(current_axes, ...
    'Box', 'On', ...
    'XGrid', 'On', ...
    'YGrid', 'On', ...
    'GridAlpha', 0.5, ...
    'FontSize', 24, ...
    'NextPlot','replacechildren', ...
    'TickLabelInterpreter','latex', ...
    'ColorOrder', colors)

current_axes.XLabel.Color = [0 0 0];
current_axes.YLabel.Color = [0 0 0];
end


