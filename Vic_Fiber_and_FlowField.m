% Plot fiber over the flow field (initial and contact positions)

clear; close all; clc;

% load flow field simulation data
load(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\Data_Give_to' ...
    '_Zhibo_20230223\input_data\Mid_plane_data.mat'])

fig = figure('color', 'w');set(gcf, 'Position', [100 100 1000 800]);

t = tiledlayout(2,1,'TileSpacing','none','Padding','tight');
nexttile

% % plot contour 
% contourf(XX*scale,YY*scale,U_norm/U_max,100,'LineStyle','none'); hold on;
% shading interp
axis equal

xlim([800 1300]*1E-6)
ylim([300 500]*1E-6)

% colormap('jet')

% % Plot streanlines
hold on
startX = 800*1E-6 * ones(1, 75);
startY = linspace(200*1E-6, 600*1E-6,75);
lineobj = streamline(XX*scale,YY*scale, UX_midplane/U_max,UY_midplane/U_max, startX, startY);
% quiver(XX*scale,YY*scale, UX_midplane/U_max,UY_midplane/U_max)
for foo = 1:length(lineobj)
    lineobj(foo).LineWidth = 0.1;
    lineobj(foo).Color = [.6, .6, .6];
end
axis off

% % Plot separatrix
hold on 
lineobj_sep = streamline(XX*scale,YY*scale, UX_midplane/U_max,UY_midplane/U_max, 800*1E-6, 399.4*1E-6);
lineobj_sep.Color = 'red';
lineobj_sep.LineWidth = 1.5;

% % Plot box edge
hold on;
rectangle('Position',[800 300 500 200]*1E-6, 'LineStyle','--')

% The position of the obstacle
obs = readVTK(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared' ...
    '\Data_Give_to_Zhibo_20230223\input_data\obstacle_beads.vtk']);
obs_beads_size = 0.6e-6; % notice here (radius)
obs_2d = obs.points(1:208, 1:2); obs_2d_center = mean(obs_2d, 1);  % 208 because there are multiple layers.
% sort the coordinates clockwise
[theta, ~] = cart2pol(obs_2d(:,1)-obs_2d_center(1), obs_2d(:,2)-obs_2d_center(2));
obs_2d = sortrows([obs_2d, theta], 3); obs_2d = obs_2d(:, 1:2);
% obstacle outer edge
obs_2d(:, 1) = (obs_2d(:, 1) - obs_2d_center(1)) * ...
    (2.5e-5 + obs_beads_size) / 2.5e-5 + obs_2d_center(1);
obs_2d(:, 2) = (obs_2d(:, 2) - obs_2d_center(2)) * ...
    (2.5e-5 + obs_beads_size) / 2.5e-5 + obs_2d_center(2);  % 2.5e-5 is the 1/3 the height of the equilateral triangle
obs_2d = [obs_2d; obs_2d(1, :)];
hold on; plot(obs_2d(:,1), obs_2d(:,2), 'k', 'LineWidth', 2);

% load fiber and plot it
% theta0: -5   y0: 0.35
snapshot = readVTK(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\simulations\theta0_m5\L_1\y0_0o35\output_data\fibers_000045.vtk']);
fiber_XY = snapshot.points(:, 1:2);
hold on
plot(fiber_XY(:,1), fiber_XY(:,2), 'Color', 'k', 'LineWidth', 3)
snapshot = readVTK(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\simulations\theta0_m5\L_1\y0_0o35\output_data\fibers_000070.vtk']);
fiber_XY = snapshot.points(:, 1:2);
hold on
plot(fiber_XY(:,1), fiber_XY(:,2), 'Color', 'k', 'LineWidth', 3)
snapshot = readVTK(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\simulations\theta0_m5\L_1\y0_0o35\output_data\fibers_000083.vtk']);
fiber_XY = snapshot.points(:, 1:2);
hold on
plot(fiber_XY(:,1), fiber_XY(:,2), 'Color', 'k', 'LineWidth', 3)
% theta0: -5   y0: 0.5
snapshot = readVTK(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\simulations\theta0_m5\L_1\y0_0o5\output_data\fibers_000045.vtk']);
fiber_XY = snapshot.points(:, 1:2);
hold on
plot(fiber_XY(:,1), fiber_XY(:,2), 'Color', 'm', 'LineWidth', 3)
snapshot = readVTK(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\simulations\theta0_m5\L_1\y0_0o5\output_data\fibers_000070.vtk']);
fiber_XY = snapshot.points(:, 1:2);
hold on
plot(fiber_XY(:,1), fiber_XY(:,2), 'Color', 'm', 'LineWidth', 3)
snapshot = readVTK(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\simulations\theta0_m5\L_1\y0_0o5\output_data\fibers_000089.vtk']);
fiber_XY = snapshot.points(:, 1:2);
hold on
plot(fiber_XY(:,1), fiber_XY(:,2), 'Color', 'm', 'LineWidth', 3)

set_plot(gcf, gca)

nexttile

% % plot contour 
% contourf(XX*scale,YY*scale,U_norm/U_max,100,'LineStyle','none'); hold on;
% shading interp
axis equal

xlim([800 1300]*1E-6)
ylim([300 500]*1E-6)

% colormap('jet')

% % Plot streanlines
hold on
startX = 800*1E-6 * ones(1, 75);
startY = linspace(200*1E-6, 600*1E-6,75);
lineobj = streamline(XX*scale,YY*scale, UX_midplane/U_max,UY_midplane/U_max, startX, startY);
% quiver(XX*scale,YY*scale, UX_midplane/U_max,UY_midplane/U_max)
for foo = 1:length(lineobj)
    lineobj(foo).LineWidth = 0.1;
    lineobj(foo).Color = [.6, .6, .6];
end
axis off

% % Plot separatrix
hold on 
lineobj_sep = streamline(XX*scale,YY*scale, UX_midplane/U_max,UY_midplane/U_max, 800*1E-6, 399.4*1E-6);
lineobj_sep.Color = 'red';
lineobj_sep.LineWidth = 1.5;

% % Plot box edge
hold on;
rectangle('Position',[800 300 500 200]*1E-6, 'LineStyle','--')

% The position of the obstacle
obs = readVTK(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared' ...
    '\Data_Give_to_Zhibo_20230223\input_data\obstacle_beads.vtk']);
obs_beads_size = 0.6e-6; % notice here (radius)
obs_2d = obs.points(1:208, 1:2); obs_2d_center = mean(obs_2d, 1);  % 208 because there are multiple layers.
% sort the coordinates clockwise
[theta, ~] = cart2pol(obs_2d(:,1)-obs_2d_center(1), obs_2d(:,2)-obs_2d_center(2));
obs_2d = sortrows([obs_2d, theta], 3); obs_2d = obs_2d(:, 1:2);
% obstacle outer edge
obs_2d(:, 1) = (obs_2d(:, 1) - obs_2d_center(1)) * ...
    (2.5e-5 + obs_beads_size) / 2.5e-5 + obs_2d_center(1);
obs_2d(:, 2) = (obs_2d(:, 2) - obs_2d_center(2)) * ...
    (2.5e-5 + obs_beads_size) / 2.5e-5 + obs_2d_center(2);  % 2.5e-5 is the 1/3 the height of the equilateral triangle
obs_2d = [obs_2d; obs_2d(1, :)];
hold on; plot(obs_2d(:,1), obs_2d(:,2), 'k', 'LineWidth', 2);

% load fiber and plot it
% theta0: -5   y0: 0.35
snapshot = readVTK(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\simulations\theta0_m5\L_1\y0_0o35\output_data\fibers_000045.vtk']);
fiber_XY = snapshot.points(:, 1:2);
hold on
plot(fiber_XY(:,1), fiber_XY(:,2), 'Color', 'k', 'LineWidth', 3)
snapshot = readVTK(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\simulations\theta0_m5\L_1\y0_0o35\output_data\fibers_000070.vtk']);
fiber_XY = snapshot.points(:, 1:2);
hold on
plot(fiber_XY(:,1), fiber_XY(:,2), 'Color', 'k', 'LineWidth', 3)
snapshot = readVTK(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\simulations\theta0_m5\L_1\y0_0o35\output_data\fibers_000083.vtk']);
fiber_XY = snapshot.points(:, 1:2);
hold on
plot(fiber_XY(:,1), fiber_XY(:,2), 'Color', 'k', 'LineWidth', 3)
% theta0: 2.5   y0: 0.35
snapshot = readVTK(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\simulations\theta0_2o5\L_1\y0_0o35\output_data\fibers_000045.vtk']);
fiber_XY = snapshot.points(:, 1:2);
hold on
plot(fiber_XY(:,1), fiber_XY(:,2), 'Color', 'g', 'LineWidth', 3)
snapshot = readVTK(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\simulations\theta0_2o5\L_1\y0_0o35\output_data\fibers_000070.vtk']);
fiber_XY = snapshot.points(:, 1:2);
hold on
plot(fiber_XY(:,1), fiber_XY(:,2), 'Color', 'g', 'LineWidth', 3)
snapshot = readVTK(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\simulations\theta0_2o5\L_1\y0_0o35\output_data\fibers_000086.vtk']);
fiber_XY = snapshot.points(:, 1:2);
hold on
plot(fiber_XY(:,1), fiber_XY(:,2), 'Color', 'g', 'LineWidth', 3)



% cb = colorbar;
% cb.Layout.Tile = 'east';
% cb.Label.String = '$U/U_{\rm max}$';
% cb.Label.Interpreter = 'LaTeX';
% cb.TickLabelInterpreter = 'LaTeX';
% cb.FontSize = 24;

set_plot(gcf, gca)

f=gcf;
set(f,'renderer','Painters');
print('-depsc2','-tiff','-r300','-painters',['F:\Processing & Results\' ...
    'FSI - Rigid Fiber &  Individual Obstacle\Figures\Paper\' ...
    'contact information vs initial condition\fiber-over-field_without_flowfield.eps'])
