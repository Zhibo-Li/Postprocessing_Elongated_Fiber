%% This is applied for the figures of Poincare maps of the fibers from Clement's simulations
% The points are plotted as 'Line' instead of 'Scatter' in the original figure.

clear; close all; clc;

plot_angle = 45;
fig = openfig(['Z:\Processing & Results\Actin Filaments in Porous Media\Figures\' ...
    'Poincare plots\Simu_poincare_maps_alpha_45.fig']);

% change the marker styles, and delete some lines
h = findobj(fig, 'Type', 'Line');
PointNumberOnLine = zeros(length(h), 1);
for i = 1:length(h)
    PointNumberOnLine(i) = numel(h(i).XData);
end

tracerSymbols = find(PointNumberOnLine <= 4);
fiberSymbols = find(PointNumberOnLine > 4 & PointNumberOnLine ~= 100);

randomNum = randperm(numel(tracerSymbols), round(0.5*numel(tracerSymbols)));
for i = tracerSymbols(1):tracerSymbols(end)
    if ismember(i, randomNum)
        delete(h(i)); continue;
    end
    h(i).XData = h(i).XData(1); h(i).YData = h(i).YData(1);
    h(i).LineStyle = "none";
    h(i).Marker = 'x';
    h(i).MarkerSize = 10;
    h(i).LineWidth = 2;

    h(i).Color = 'r'; % for 0 and 45 degree

    % if h(i).XData > h(i).YData
    %     h(i).Color = 'b';
    % else
    %     h(i).Color = 'r';
    % end
end

for i = fiberSymbols(1):fiberSymbols(end)
    % randomly delete 90% of the points
    randomNum = randperm(numel(h(i).XData), round(0.9*numel(h(i).XData)));
    h(i).XData(randomNum) = []; h(i).YData(randomNum) = [];
    h(i).LineWidth = 1;
    h(i).Marker = 'o';
    h(i).MarkerFaceColor = 'none';
    h(i).MarkerSize = 15;
end


% Experiment data for fibers
load('Z:\Processing & Results\Actin Filaments in Porous Media\Figures\Poincare plots\Poincare_Map_data_202306.mat');

Obj_Mag = 0.1; % um/pixel

theFlAng = Info.FlowAngle;
[C, ia, ic] = unique(theFlAng,'stable');

for ii = 1: length(C)

    Lattice_in_all = []; Lattice_out_all = []; L_toPlot_all = [];
    current_deg = C(ii);

    if current_deg == plot_angle
       
    % the contour lengths
    L_all = Info.L * Obj_Mag;

    % in & out position in a unit cell
    current_deg_index = find(theFlAng == current_deg);
    for jj = 1:2:length(current_deg_index)
        Lattice_in = Info.map{1, current_deg_index(jj)};

        Lattice_out = [[Lattice_in(2:end, 1);nan], Lattice_in(:, 2:3)];
        out_in_diff = Lattice_out(:, 2) - Lattice_in(:, 2);

        Lattice_in = Lattice_in(out_in_diff==0)';
        Lattice_out = Lattice_out(out_in_diff==0)';

        Lattice_in_all = [Lattice_in_all, Lattice_in];
        Lattice_out_all = [Lattice_out_all, Lattice_out];
        L_toPlot_all = [L_toPlot_all, L_all(current_deg_index(jj)) * ones(1, numel(Lattice_in))];

    end

    randomNum = randperm(numel(Lattice_in_all), round(0.6*numel(Lattice_in_all)));
    Lattice_in_all(randomNum) = []; Lattice_out_all(randomNum) = []; L_toPlot_all(randomNum) = [];
    scatter(Lattice_in_all, Lattice_out_all, 100, L_toPlot_all, ...
        'MarkerFaceColor','none','Marker','^','LineWidth', 1);
    cmocean('amp'); clim([0 50]);

    else
        continue
    end
end

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['Z:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Poincare plots\Poincare_', ...
    num2str(plot_angle),'def_Exp_Simu_Tracer_groupTracers.eps']);



%% 2025
% This is applied for the figures of Poincare maps of the fibers from Clement's simulations
% The points are plotted as 'Line' instead of 'Scatter' in the original figure.
% Simu. and exp. are separately plotted

%% Simulation replot
clear; close all; clc;

plot_angle = 0;
fig = openfig(['Z:\Processing & Results\Actin Filaments in Porous Media\Figures\' ...
    'Poincare plots\Simu_poincare_maps_alpha_0.fig']);

% change the marker styles, and delete some lines
h = findobj(fig, 'Type', 'Line');
PointNumberOnLine = zeros(length(h), 1);
for i = 1:length(h)
    PointNumberOnLine(i) = numel(h(i).XData);
end

tracerSymbols = find(PointNumberOnLine <= 4);
fiberSymbols = find(PointNumberOnLine > 4 & PointNumberOnLine ~= 100);

% randomNum = randperm(numel(tracerSymbols), round(0.5*numel(tracerSymbols)));
randomNum = [];
for i = tracerSymbols(1):tracerSymbols(end)
    if ismember(i, randomNum)
        delete(h(i)); continue;
    end
    h(i).XData = h(i).XData(1); h(i).YData = h(i).YData(1);
    h(i).LineStyle = "none";
    h(i).Marker = 'x';
    h(i).MarkerSize = 10;
    h(i).LineWidth = 2;

    h(i).Color = 'b'; % for 0 and 45 degree

    % if h(i).XData > h(i).YData
    %     h(i).Color = 'b';
    % else
    %     h(i).Color = 'r';
    % end
end

for i = fiberSymbols(1):fiberSymbols(end)
    % randomly delete 80% of the points (for 0 degrees)
    randomNum = randperm(numel(h(i).XData), round(0.8*numel(h(i).XData)));
    % randomNum = [];
    h(i).XData(randomNum) = []; h(i).YData(randomNum) = [];
    h(i).LineWidth = 1;
    h(i).Marker = 'o';
    h(i).MarkerFaceColor = 'none';
    h(i).MarkerSize = 15;
end

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['Z:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Poincare plots\2025 New\Poincare_', ...
    num2str(plot_angle),'deg_Simu_groupTracers.eps']);


%% Experiment replot
clear; close all; clc;

plot_angle = 45;
fig = openfig(['Z:\Processing & Results\Actin Filaments in Porous Media\Figures\' ...
    'Poincare plots\Simu_poincare_maps_alpha_45.fig']);

% change the marker styles, and delete some lines
h = findobj(fig, 'Type', 'Line');
PointNumberOnLine = zeros(length(h), 1);
for i = 1:length(h)
    PointNumberOnLine(i) = numel(h(i).XData);
end

tracerSymbols = find(PointNumberOnLine <= 4);
fiberSymbols = find(PointNumberOnLine > 4 & PointNumberOnLine ~= 100);

% randomNum = randperm(numel(tracerSymbols), round(0.5*numel(tracerSymbols)));
randomNum = [];
for i = tracerSymbols(1):tracerSymbols(end)
    if ismember(i, randomNum)
        delete(h(i)); continue;
    end
    h(i).XData = h(i).XData(1); h(i).YData = h(i).YData(1);
    h(i).LineStyle = "none";
    h(i).Marker = 'x';
    h(i).MarkerSize = 10;
    h(i).LineWidth = 2;

    h(i).Color = 'r'; % for 0 and 45 degree

    % if h(i).XData > h(i).YData
    %     h(i).Color = 'b';
    % else
    %     h(i).Color = 'r';
    % end
end

for i = fiberSymbols(1):fiberSymbols(end)
    % randomly delete 90% of the points
    % randomNum = randperm(numel(h(i).XData), round(0.9*numel(h(i).XData)));
    % h(i).XData(randomNum) = []; h(i).YData(randomNum) = [];
    h(i).XData = []; h(i).YData = [];
    h(i).LineWidth = 1;
    h(i).Marker = 'o';
    h(i).MarkerFaceColor = 'none';
    h(i).MarkerSize = 15;
end

load('Z:\Processing & Results\Actin Filaments in Porous Media\Figures\Poincare plots\Poincare_Map_data_202306.mat');

Obj_Mag = 0.1; % um/pixel

theFlAng = Info.FlowAngle;
[C, ia, ic] = unique(theFlAng,'stable');

for ii = 1: length(C)

    Lattice_in_all = []; Lattice_out_all = []; L_toPlot_all = [];
    current_deg = C(ii);

    if current_deg == plot_angle
       
    % the contour lengths
    L_all = Info.L * Obj_Mag;

    % in & out position in a unit cell
    current_deg_index = find(theFlAng == current_deg);
    for jj = 1:2:length(current_deg_index)
        Lattice_in = Info.map{1, current_deg_index(jj)};

        Lattice_out = [[Lattice_in(2:end, 1);nan], Lattice_in(:, 2:3)];
        out_in_diff = Lattice_out(:, 2) - Lattice_in(:, 2);

        Lattice_in = Lattice_in(out_in_diff==0)';
        Lattice_out = Lattice_out(out_in_diff==0)';

        Lattice_in_all = [Lattice_in_all, Lattice_in];
        Lattice_out_all = [Lattice_out_all, Lattice_out];
        L_toPlot_all = [L_toPlot_all, L_all(current_deg_index(jj)) * ones(1, numel(Lattice_in))];

    end

    % randomNum = randperm(numel(Lattice_in_all), round(0.6*numel(Lattice_in_all)));
    randomNum = [];
    Lattice_in_all(randomNum) = []; Lattice_out_all(randomNum) = []; L_toPlot_all(randomNum) = [];
    scatter(Lattice_in_all, Lattice_out_all, 150, L_toPlot_all, ...
        'MarkerFaceColor','none','Marker','^','LineWidth', 1);
    cmocean('amp'); clim([0 50]);

    else
        continue
    end
end

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['Z:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Poincare plots\2025 New\Poincare_', ...
    num2str(plot_angle),'deg_Exp_groupTracers.eps']);









%% Try different colorbar for fiber lengths
%% Experiment replot (different colorbar)
clear; close all; clc;

Vic_cmap = slanCM(169); % slanCM169: neon
% reverse the colormap
% Vic_cmap = Vic_cmap(end:-1:1, :);

plot_angle = 45;
fig = openfig(['Z:\Processing & Results\Actin Filaments in Porous Media\Figures\' ...
    'Poincare plots\Simu_poincare_maps_alpha_45.fig']);

% change the marker styles, and delete some lines
h = findobj(fig, 'Type', 'Line');
PointNumberOnLine = zeros(length(h), 1);
for i = 1:length(h)
    PointNumberOnLine(i) = numel(h(i).XData);
end

tracerSymbols = find(PointNumberOnLine <= 4);
fiberSymbols = find(PointNumberOnLine > 4 & PointNumberOnLine ~= 100);

% randomNum = randperm(numel(tracerSymbols), round(0.5*numel(tracerSymbols)));
randomNum = [];
for i = tracerSymbols(1):tracerSymbols(end)
    if ismember(i, randomNum)
        delete(h(i)); continue;
    end
    h(i).XData = h(i).XData(1); h(i).YData = h(i).YData(1);
    h(i).LineStyle = "none";
    h(i).Marker = 'x';
    h(i).MarkerSize = 10;
    h(i).LineWidth = 2;

    h(i).Color = [0.2 0.2 0.2];
    % h(i).Color = 'k';

end

for i = fiberSymbols(1):fiberSymbols(end)
    % randomly delete 90% of the points
    % randomNum = randperm(numel(h(i).XData), round(0.9*numel(h(i).XData)));
    % h(i).XData(randomNum) = []; h(i).YData(randomNum) = [];
    h(i).XData = []; h(i).YData = [];
    h(i).LineWidth = 1;
    h(i).Marker = 'o';
    h(i).MarkerFaceColor = 'none';
    h(i).MarkerSize = 15;
end

load('Z:\Processing & Results\Actin Filaments in Porous Media\Figures\Poincare plots\Poincare_Map_data_202306.mat');

Obj_Mag = 0.1; % um/pixel

theFlAng = Info.FlowAngle;
[C, ia, ic] = unique(theFlAng,'stable');

for ii = 1: length(C)

    Lattice_in_all = []; Lattice_out_all = []; L_toPlot_all = [];
    current_deg = C(ii);

    if current_deg == plot_angle
       
    % the contour lengths
    L_all = Info.L * Obj_Mag;

    % in & out position in a unit cell
    current_deg_index = find(theFlAng == current_deg);
    for jj = 1:2:length(current_deg_index)
        Lattice_in = Info.map{1, current_deg_index(jj)};

        Lattice_out = [[Lattice_in(2:end, 1);nan], Lattice_in(:, 2:3)];
        out_in_diff = Lattice_out(:, 2) - Lattice_in(:, 2);

        Lattice_in = Lattice_in(out_in_diff==0)';
        Lattice_out = Lattice_out(out_in_diff==0)';

        Lattice_in_all = [Lattice_in_all, Lattice_in];
        Lattice_out_all = [Lattice_out_all, Lattice_out];
        L_toPlot_all = [L_toPlot_all, L_all(current_deg_index(jj)) * ones(1, numel(Lattice_in))];

    end

    % randomNum = randperm(numel(Lattice_in_all), round(0.6*numel(Lattice_in_all)));
    randomNum = [];
    Lattice_in_all(randomNum) = []; Lattice_out_all(randomNum) = []; L_toPlot_all(randomNum) = [];
    scatter(Lattice_in_all, Lattice_out_all, 120, L_toPlot_all, ...
        'MarkerFaceColor','none','Marker','^','LineWidth', 1);
    colormap(Vic_cmap); clim([0 50]); 
    % cmocean('-thermal'); clim([0 50]);

    else
        continue
    end
end

% colorbar;
set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['Z:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Poincare plots\2025 New\Poincare_', ...
    num2str(plot_angle),'deg_Exp_groupTracers_diffColorbar-neon.eps']);



%% Simulation replot (different colorbar)
clear; close all; clc;

Vic_cmap = slanCM(169); % slanCM169: neon
% reverse the colormap
% Vic_cmap = Vic_cmap(end:-1:1, :);

% Get the colormap data
% cmap = cmocean('-thermal');
cmap = Vic_cmap;
max_length = 50;
min_length = 0;
shorter_length = 12;
longer_length = 48;
% Determine the color index
color_index_short = round((shorter_length - min_length) / (max_length - min_length) * (size(cmap, 1) - 1)) + 1;
color_index_long = round((longer_length - min_length) / (max_length - min_length) * (size(cmap, 1) - 1)) + 1;
% Get the color value
color_value_short = cmap(color_index_short, :);
color_value_long = cmap(color_index_long, :);

plot_angle = 30;
fig = openfig(['Z:\Processing & Results\Actin Filaments in Porous Media\Figures\' ...
    'Poincare plots\Simu_poincare_maps_alpha_30.fig']);

% change the marker styles, and delete some lines
h = findobj(fig, 'Type', 'Line');
PointNumberOnLine = zeros(length(h), 1);
for i = 1:length(h)
    PointNumberOnLine(i) = numel(h(i).XData);
end

tracerSymbols = find(PointNumberOnLine <= 4);
fiberSymbols = find(PointNumberOnLine > 4 & PointNumberOnLine ~= 100);

% randomNum = randperm(numel(tracerSymbols), round(0.5*numel(tracerSymbols)));
randomNum = [];
for i = tracerSymbols(1):tracerSymbols(end)
    if ismember(i, randomNum)
        delete(h(i)); continue;
    end
    h(i).XData = h(i).XData(1); h(i).YData = h(i).YData(1);
    h(i).LineStyle = "none";
    h(i).Marker = 'x';
    h(i).MarkerSize = 10;
    h(i).LineWidth = 2;

    h(i).Color = [0.2 0.2 0.2];
    % h(i).Color = 'k';

end

for i = fiberSymbols(1):fiberSymbols(end)
    % randomly delete 80% of the points (for 0 degrees)
    % randomNum = randperm(numel(h(i).XData), round(0.8*numel(h(i).XData)));
    randomNum = [];
    h(i).XData(randomNum) = []; h(i).YData(randomNum) = [];
    h(i).LineWidth = 1;
    h(i).Marker = 'o';
    h(i).MarkerFaceColor = 'none';
    h(i).MarkerSize = 10;
    if h(i).Color == [0.85, 0.66, 0.59] % this short fibers
        h(i).Color = color_value_short;
    else
        h(i).Color = color_value_long;
    end
end

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['Z:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Poincare plots\2025 New\Poincare_', ...
    num2str(plot_angle),'deg_Simu_groupTracers_diffColorbar-neon.eps']);



%% Plot only the colorbar
figure('color', 'w'); set(gcf, 'Position', [100 100 200 380]);

Vic_cmap = slanCM(169); % slanCM169: neon
% reverse the colormap
% Vic_cmap = Vic_cmap(end:-1:1, :);
colormap(Vic_cmap); clim([0 50])
c = colorbar;
c.Label.String = '$L\ (\mathrm{\mu m})$';
c.Label.Interpreter = 'LaTeX';
c.TickLabelInterpreter = 'LaTeX';
c.FontSize = 20;
c.Location = 'west';
axis off

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['Z:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Poincare plots\2025 New\The_colorbar_Neon.eps']);