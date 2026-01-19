%% 1st version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

dt = 0.1; % time step, unit: s

% Geometrical and mechanical properties of the fiber
R_fib = 3e-7; % fiber radius, unit: m
B = 8e-25; % Bending rigidity
S = 4*B/R_fib^2;
ks = 100*S/(2*R_fib);

u_max = 200e-6;
mu = 6.1e-3;
f_bead = 6*pi*mu*R_fib*u_max;

% Geometrical properties of the array
R_pillar = 8.3 * 1e-6; % pillar radius, unit: m
C2C_pillar = 30 * 1e-6; % m (delta_x & delta_y)

time_step = 0.1; % s

X_pillar = -30*C2C_pillar:C2C_pillar:100*C2C_pillar; Y_pillar = 30*C2C_pillar:-C2C_pillar:-100*C2C_pillar;
[X_pillar, Y_pillar] = meshgrid(X_pillar, Y_pillar);

theta = -35; % angle in degrees (clockwise)
Rotate_matrix = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; % rotation matrix

XY_pillar_rotated = [X_pillar(:), Y_pillar(:)]*Rotate_matrix;
% % figure('Position', [100, 100, 800, 800], 'Color', 'w');
% % hold on;
% % viscircles(XY_pillar_rotated, R_pillar*ones(size(XY_pillar_rotated, 1), 1), 'Color', 'r');
% % axis equal; axis off

simu_fiber_parent_path_char = 'Z:\Collaboration - LadHyX\Actin in PAs\SI videos data\fibers';
simu_fiber_parent_path = dir(simu_fiber_parent_path_char);
simu_fiber_parent_path(ismember({simu_fiber_parent_path.name}, {'.', '..'})) = [];

for folder_idx = 1:length(simu_fiber_parent_path)

    current_fiber_path = fullfile(simu_fiber_parent_path(folder_idx).folder, simu_fiber_parent_path(folder_idx).name);
    if ~isfolder(current_fiber_path)
        continue;
    end

    if contains(current_fiber_path, 'short')
        case_no = 1;
    elseif contains(current_fiber_path, 'medium')
        case_no = 2;
    elseif contains(current_fiber_path, 'long')
        case_no = 3;
    end

    fileinfo = dir(fullfile(current_fiber_path, '\*.vtk'));
    if isempty(fileinfo)
        continue;
    end
    
    if case_no == 1
        picked_snaps_short =  [0 27 57 87 118 150 185 212 242 270 300 327 354:35:1480]+1;
        XY_rotated_short = zeros(length(picked_snaps_short), length(readVTK(fullfile(fileinfo(1).folder, fileinfo(1).name)).points), 2);
        n = 1;
        for ii = picked_snaps_short
            snapshot = readVTK(fullfile(fileinfo(ii).folder, fileinfo(ii).name));
            XY_current = snapshot.points(:, 1:2);
            XY_rotated_short(n, :, :) = XY_current*Rotate_matrix;
            n = n + 1;
        end

    elseif case_no == 2
        picked_snaps_medium = [0 48 160 240 310:88:1100]+1;
        XY_rotated_medium = zeros(length(picked_snaps_medium), length(readVTK(fullfile(fileinfo(1).folder, fileinfo(1).name)).points), 2);
        n = 1;
        for ii = picked_snaps_medium
            snapshot = readVTK(fullfile(fileinfo(ii).folder, fileinfo(ii).name));
            XY_current = snapshot.points(:, 1:2);
            XY_rotated_medium(n, :, :) = XY_current*Rotate_matrix;
            n = n + 1;
        end

    elseif case_no == 3
        picked_snaps_long = [10 70 130 200:70:500]+1;
        XY_rotated_long = zeros(length(picked_snaps_long), length(readVTK(fullfile(fileinfo(1).folder, fileinfo(1).name)).points), 2);
        n = 1;
        for ii = picked_snaps_long
            snapshot = readVTK(fullfile(fileinfo(ii).folder, fileinfo(ii).name));
            XY_current = snapshot.points(:, 1:2);
            XY_rotated_long(n, :, :) = XY_current*Rotate_matrix;
            n = n + 1;
        end
    end



end

% Plot the chronophotography (short fiber)
figure('Position', [100, 100, 2000, 600], 'Color', 'w');
Fiber_num = size(XY_rotated_short, 1);
the_RGB_color_short = [255, 18, 18]/255;
the_HSV_color_short = rgb2hsv(the_RGB_color_short);
HSV_colors_short = repmat(the_HSV_color_short, Fiber_num, 1);
HSV_colors_short(:, 2) = linspace(0.4, the_HSV_color_short(2), Fiber_num)';
RGB_colors_short = hsv2rgb(HSV_colors_short);
for ii = 1:Fiber_num
    plot(XY_rotated_short(ii, :, 1), XY_rotated_short(ii, :, 2), 'Color', RGB_colors_short(ii, :), 'LineWidth', 3);
    hold on;
end
% Plot the chronophotography (medium fiber)
Fiber_num = size(XY_rotated_medium, 1);
the_RGB_color_medium = [12, 178, 34]/255;
the_HSV_color_medium = rgb2hsv(the_RGB_color_medium);
HSV_colors_medium = repmat(the_HSV_color_medium, Fiber_num, 1);
HSV_colors_medium(:, 2) = linspace(0.4, the_HSV_color_medium(2), Fiber_num)';
RGB_colors_medium = hsv2rgb(HSV_colors_medium);
for ii = 1:Fiber_num
    plot(XY_rotated_medium(ii, :, 1), XY_rotated_medium(ii, :, 2), 'Color', RGB_colors_medium(ii, :), 'LineWidth', 3);
    hold on;
end
% Plot the chronophotography (long fiber)
Fiber_num = size(XY_rotated_long, 1);
the_RGB_color_long = [14, 35, 202]/255;
the_HSV_color_long = rgb2hsv(the_RGB_color_long);
HSV_colors_long = repmat(the_HSV_color_long, Fiber_num, 1);
HSV_colors_long(:, 2) = linspace(0.4, the_HSV_color_long(2), Fiber_num)';
RGB_colors_long = hsv2rgb(HSV_colors_long);
for ii = 1:Fiber_num
    plot(XY_rotated_long(ii, :, 1), XY_rotated_long(ii, :, 2), 'Color', RGB_colors_long(ii, :), 'LineWidth', 3);
    hold on;
end


XY_pillar_rotated(XY_pillar_rotated(:, 2) < -0.5e-4, :) = [];
XY_pillar_rotated(XY_pillar_rotated(:, 2) > 1.5e-4, :) = [];
XY_pillar_rotated(XY_pillar_rotated(:, 1) < -0.5e-4, :) = [];
XY_pillar_rotated(XY_pillar_rotated(:, 1) > 6.5e-4, :) = [];
viscircles(XY_pillar_rotated, R_pillar*ones(size(XY_pillar_rotated, 1), 1), 'Color', 'k', 'LineWidth', 0.3);
axis equal; 
xlim([-0.3 4.9]*1e-4); ylim([-0.2 1]*1e-4);
box on;
set(gca, 'xtick', [], 'ytick', []);
set(gca,'LineWidth',1);
set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',[simu_fiber_parent_path_char, filesep, 'Three_fibers_chronophotography.eps']);










%% 2nd version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

dt = 0.1; % time step, unit: s

% Geometrical and mechanical properties of the fiber
R_fib = 3e-7; % fiber radius, unit: m
B = 8e-25; % Bending rigidity
S = 4*B/R_fib^2;
ks = 100*S/(2*R_fib);

u_max = 200e-6;
mu = 6.1e-3;
f_bead = 6*pi*mu*R_fib*u_max;

% Geometrical properties of the array
R_pillar = 8.3 * 1e-6; % pillar radius, unit: m
C2C_pillar = 30 * 1e-6; % m (delta_x & delta_y)

time_step = 0.1; % s

X_pillar = -30*C2C_pillar:C2C_pillar:100*C2C_pillar; Y_pillar = 30*C2C_pillar:-C2C_pillar:-100*C2C_pillar;
[X_pillar, Y_pillar] = meshgrid(X_pillar, Y_pillar);

theta = -35; % angle in degrees (clockwise)
Rotate_matrix = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; % rotation matrix

XY_pillar_rotated = [X_pillar(:), Y_pillar(:)]*Rotate_matrix;
% % figure('Position', [100, 100, 800, 800], 'Color', 'w');
% % hold on;
% % viscircles(XY_pillar_rotated, R_pillar*ones(size(XY_pillar_rotated, 1), 1), 'Color', 'r');
% % axis equal; axis off

simu_fiber_parent_path_char = 'Z:\Collaboration - LadHyX\Actin in PAs\SI videos data\fibers';
simu_fiber_parent_path = dir(simu_fiber_parent_path_char);
simu_fiber_parent_path(ismember({simu_fiber_parent_path.name}, {'.', '..'})) = [];

for folder_idx = 1:length(simu_fiber_parent_path)

    current_fiber_path = fullfile(simu_fiber_parent_path(folder_idx).folder, simu_fiber_parent_path(folder_idx).name);
    if ~isfolder(current_fiber_path)
        continue;
    end

    if contains(current_fiber_path, 'short')
        case_no = 1;
    elseif contains(current_fiber_path, 'medium')
        case_no = 2;
    elseif contains(current_fiber_path, 'long')
        case_no = 3;
    end

    fileinfo = dir(fullfile(current_fiber_path, '\*.vtk'));
    if isempty(fileinfo)
        continue;
    end
    
    if case_no == 1
        picked_snaps_short =  [0 27 57 87 118 150 185 212 242 270 300 327 357:36:1480]+1;
        XY_rotated_short = zeros(length(picked_snaps_short), length(readVTK(fullfile(fileinfo(1).folder, fileinfo(1).name)).points), 2);
        n = 1;
        for ii = picked_snaps_short
            snapshot = readVTK(fullfile(fileinfo(ii).folder, fileinfo(ii).name));
            XY_current = snapshot.points(:, 1:2);
            XY_rotated_short(n, :, :) = XY_current*Rotate_matrix;
            n = n + 1;
        end

    elseif case_no == 2
        picked_snaps_medium = [0 48 160 240 327:89:1100]+1;
        XY_rotated_medium = zeros(length(picked_snaps_medium), length(readVTK(fullfile(fileinfo(1).folder, fileinfo(1).name)).points), 2);
        n = 1;
        for ii = picked_snaps_medium
            snapshot = readVTK(fullfile(fileinfo(ii).folder, fileinfo(ii).name));
            XY_current = snapshot.points(:, 1:2);
            XY_rotated_medium(n, :, :) = XY_current*Rotate_matrix;
            n = n + 1;
        end

    elseif case_no == 3
        picked_snaps_long = [10 70 130 219:83:510]+1;
        XY_rotated_long = zeros(length(picked_snaps_long), length(readVTK(fullfile(fileinfo(1).folder, fileinfo(1).name)).points), 2);
        n = 1;
        for ii = picked_snaps_long
            snapshot = readVTK(fullfile(fileinfo(ii).folder, fileinfo(ii).name));
            XY_current = snapshot.points(:, 1:2);
            XY_rotated_long(n, :, :) = XY_current*Rotate_matrix;
            n = n + 1;
        end
    end



end

% Plot the chronophotography (short fiber)
figure('Position', [100, 100, 2000, 600], 'Color', 'w');
Fiber_num = size(XY_rotated_short, 1);
the_RGB_color_short = [255, 18, 18]/255;
the_HSV_color_short = rgb2hsv(the_RGB_color_short);
HSV_colors_short = repmat(the_HSV_color_short, Fiber_num, 1);
HSV_colors_short(:, 2) = linspace(0.4, the_HSV_color_short(2), Fiber_num)';
RGB_colors_short = hsv2rgb(HSV_colors_short);
for ii = 1:Fiber_num
    plot(XY_rotated_short(ii, :, 1), XY_rotated_short(ii, :, 2), 'Color', RGB_colors_short(ii, :), 'LineWidth', 2);
    hold on;
end
XY_pillar_rotated(XY_pillar_rotated(:, 2) < -0.5e-4, :) = [];
XY_pillar_rotated(XY_pillar_rotated(:, 2) > 1.5e-4, :) = [];
XY_pillar_rotated(XY_pillar_rotated(:, 1) < -0.5e-4, :) = [];
XY_pillar_rotated(XY_pillar_rotated(:, 1) > 6.5e-4, :) = [];
viscircles(XY_pillar_rotated, R_pillar*ones(size(XY_pillar_rotated, 1), 1), 'Color', 'k', 'LineWidth', 0.3);
axis equal; 
xlim([-0.1 5.1]*1e-4); ylim([-0.1 0.7]*1e-4);
box on; hold off
set(gca, 'xtick', [], 'ytick', []);
set(gca,'LineWidth',1);
set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',[simu_fiber_parent_path_char, filesep, 'Short_fibers_chronophotography_2.eps']);



% % Plot the chronophotography (medium fiber)
% figure('Position', [100, 100, 2000, 600], 'Color', 'w');
% Fiber_num = size(XY_rotated_medium, 1);
% the_RGB_color_medium = [12, 178, 34]/255;
% the_HSV_color_medium = rgb2hsv(the_RGB_color_medium);
% HSV_colors_medium = repmat(the_HSV_color_medium, Fiber_num, 1);
% HSV_colors_medium(:, 2) = linspace(0.4, the_HSV_color_medium(2), Fiber_num)';
% RGB_colors_medium = hsv2rgb(HSV_colors_medium);
% for ii = 1:Fiber_num
%     plot(XY_rotated_medium(ii, :, 1), XY_rotated_medium(ii, :, 2), 'Color', RGB_colors_medium(ii, :), 'LineWidth', 2);
%     hold on;
% end
% XY_pillar_rotated(XY_pillar_rotated(:, 2) < -0.5e-4, :) = [];
% XY_pillar_rotated(XY_pillar_rotated(:, 2) > 1.5e-4, :) = [];
% XY_pillar_rotated(XY_pillar_rotated(:, 1) < -0.5e-4, :) = [];
% XY_pillar_rotated(XY_pillar_rotated(:, 1) > 6.5e-4, :) = [];
% viscircles(XY_pillar_rotated, R_pillar*ones(size(XY_pillar_rotated, 1), 1), 'Color', 'k', 'LineWidth', 0.3);
% axis equal; 
% xlim([-0.3 4.9]*1e-4); ylim([0 1]*1e-4);
% box on; hold off
% set(gca, 'xtick', [], 'ytick', []);
% set(gca,'LineWidth',1);
% set(gcf,'renderer','Painters');
% print('-depsc2','-tiff','-r100','-vector',[simu_fiber_parent_path_char, filesep, 'Medium_fibers_chronophotography.eps']);



% Plot the chronophotography (long fiber)
figure('Position', [100, 100, 2000, 600], 'Color', 'w');
Fiber_num = size(XY_rotated_long, 1);
the_RGB_color_long = [14, 35, 202]/255;
the_HSV_color_long = rgb2hsv(the_RGB_color_long);
HSV_colors_long = repmat(the_HSV_color_long, Fiber_num, 1);
HSV_colors_long(:, 2) = linspace(0.4, the_HSV_color_long(2), Fiber_num)';
RGB_colors_long = hsv2rgb(HSV_colors_long);
for ii = 1:Fiber_num
    plot(XY_rotated_long(ii, :, 1), XY_rotated_long(ii, :, 2), 'Color', RGB_colors_long(ii, :), 'LineWidth', 2);
    hold on;
end
XY_pillar_rotated(XY_pillar_rotated(:, 2) < -0.5e-4, :) = [];
XY_pillar_rotated(XY_pillar_rotated(:, 2) > 1.5e-4, :) = [];
XY_pillar_rotated(XY_pillar_rotated(:, 1) < -0.5e-4, :) = [];
XY_pillar_rotated(XY_pillar_rotated(:, 1) > 6.5e-4, :) = [];
viscircles(XY_pillar_rotated, R_pillar*ones(size(XY_pillar_rotated, 1), 1), 'Color', 'k', 'LineWidth', 0.3);
axis equal; 
xlim([-0.1 5.1]*1e-4); ylim([-0.2 0.6]*1e-4);
box on; hold off
set(gca, 'xtick', [], 'ytick', []);
set(gca,'LineWidth',1);
set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',[simu_fiber_parent_path_char, filesep, 'Long_fibers_chronophotography_2.eps']);



