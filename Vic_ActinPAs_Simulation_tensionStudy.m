%% Calculate the tension in the fiber based on the simulation results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

% Geometrical and mechanical properties of the fiber
R_fib = 3e-7; % fiber radius, unit: m
B = 8e-25; % Bending rigidity
S = 4*B/R_fib^2;
ks = 100*S/(2*R_fib);

% Geometrical properties of the array
R_pillar = 8.3 * 1e-6; % pillar radius, unit: m
C2C_pillar = 30 * 1e-6; % m (delta_x & delta_y)

time_step = 0.1; % s

data_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Simulation\Simulations_LubricationModel'];

parent_path = 'E:\Experimental Data (RAW)\Simulation\data';
sub1_path = dir(parent_path);
% only keep the folders whose name start with "n*"
sub1_path = sub1_path(arrayfun(@(x) startsWith(x.name, 'n'), sub1_path));
for sub1Path_i = 1:length(sub1_path)
    current_beadsNum = sub1_path(sub1Path_i).name;
    beads_num = str2double(current_beadsNum(2:end));

    L_0 = 0.6 * beads_num * 1e-6;

    sub2_path = dir(fullfile(parent_path, current_beadsNum));

    XY = cell(1500, length(sub2_path)-5); % 1500 is the maximum number of snapshots (roughly)
    f_tension_vector = cell(1500, length(sub2_path)-5); % 1500 is the maximum number of snapshots (roughly)
    f_tension_magnitude = cell(1500, length(sub2_path)-5); % 1500 is the maximum number of snapshots (roughly)
    ContactBeads_index = cell(1500, length(sub2_path)-5); % 1500 is the maximum number of snapshots (roughly)

    % only keep the folders whose name start with "sim*"
    sub2_path = sub2_path(arrayfun(@(x) startsWith(x.name, 'sim'), sub2_path));
    for sub2Path_i = 1:length(sub2_path)
        current_sim = sub2_path(sub2Path_i).name;
        sim_no = str2double(current_sim(11:end));

        fileinfo = dir(fullfile(parent_path, current_beadsNum, current_sim, 'output_data\*.vtk'));

        for ii = 1:length(fileinfo)
            snapshot = readVTK(fullfile(fileinfo(ii).folder, fileinfo(ii).name));
            XY_current = snapshot.points(:, 1:2);

            XY_current_in_lattice = mod(XY_current-C2C_pillar/2, C2C_pillar);
            dist2origin = sqrt((XY_current_in_lattice(:,1)-C2C_pillar/2).^2 + (XY_current_in_lattice(:,2)-C2C_pillar/2).^2);
            R_ref_Clement = R_pillar + 1.1*R_fib; % Clement's
            R_ref_Zhibo = R_pillar + R_fib; % Zhibo's
            ContactBeads_index_current = dist2origin < R_ref_Clement;

            ContactBeads_index{ii, sub2Path_i} = ContactBeads_index_current;

            % % plot XY_current_in_lattice in the unit cell with the pillars, each dot represents a bead with a radius of R_fib
            % figure('Position', [100, 100, 800, 800], 'Color', 'w');
            % for foo = 1:size(XY_current_in_lattice, 1)
            %     if ContactBeads_index_current(foo)
            %         viscircles(XY_current_in_lattice(foo, :)/C2C_pillar, R_fib/C2C_pillar, 'Color', 'm', 'LineStyle', '-'); hold on;
            %     else
            %         viscircles(XY_current_in_lattice(foo, :)/C2C_pillar, R_fib/C2C_pillar, 'Color', 'b', 'LineStyle', '-'); hold on;
            %     end
            % end
            % viscircles(XY_current_in_lattice(1, :)/C2C_pillar, R_fib/C2C_pillar, 'Color', 'r', 'LineStyle', '-'); hold on;
            % viscircles([0.5 0.5], R_pillar/C2C_pillar, 'Color', 'r', 'LineStyle', '--');
            % axis equal; axis off
            % xlim([0 1]); ylim([0 1]);
            
            XY{ii, sub2Path_i} = XY_current;

            vec_im = diff(XY_current);
            l_im = vecnorm(vec_im')';
            vec_t = vec_im./l_im;
            f_tension_vector{ii, sub2Path_i} = -ks*(l_im - 2*R_fib).*vec_t;
            f_tension_magnitude{ii, sub2Path_i} = ks*(l_im - 2*R_fib); % Notice the sign

        end
    end

    % remove the empty cell
    XY = XY(~cellfun('isempty', XY));
    ContactBeads_index = ContactBeads_index(~cellfun('isempty', ContactBeads_index));
    f_tension_vector = f_tension_vector(~cellfun('isempty', f_tension_vector));
    f_tension_magnitude = f_tension_magnitude(~cellfun('isempty', f_tension_magnitude));

    if ~exist("data_mother_save_path", 'dir')
        mkdir(data_mother_save_path);
    end

    save([data_mother_save_path, filesep, 'FiberLength=', num2str(L_0*1e6), 'um_TensionStudy.mat'], ...
        'XY', 'f_tension_vector', 'f_tension_magnitude', 'L_0', 'ContactBeads_index');

end





%% Load saved tension data and plot the results (in unit cell)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

% Geometrical properties of the array
R_pillar = 8.3 * 1e-6; % pillar radius, unit: m
C2C_pillar = 30 * 1e-6; % m (delta_x & delta_y)
delta_t = 0.1; % s

R_fib = 3e-7; % fiber radius, unit: m
u_max = 200e-6;
mu = 6.1e-3;
f_bead = 6*pi*mu*R_fib*u_max;

Xedges = 0:0.02:1; Yedges = 0:0.02:1; % for the histogram

data_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Simulation\Simulations_LubricationModel'];
saved_data = dir(data_mother_save_path);
% only keep the .mat files inluding the '_TensionStudy'
saved_data = saved_data(arrayfun(@(x) contains(x.name, '_TensionStudy.mat'), saved_data)); 

figure_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Simulation\Tension analysis'];

L_0_group = zeros(1, length(saved_data));
f_tension_NonContact_normalized_magnitude_in_grid_mean_group = zeros(length(Xedges)-1, length(Yedges)-1, length(saved_data));
for ii = 1:length(saved_data)

    load(fullfile(saved_data(ii).folder, saved_data(ii).name));

    L_0 = extractBetween(saved_data(ii).name, 'FiberLength=', 'um_TensionStudy.mat'); % unit: um
    bead_num = size(XY{1, 1}, 1);

    % Apply cellfun to shift contact beads index to align with the tension data
    % Exclude the tension calculated for beads that are in contact with the pillars
    ContactBeads_index_aligned = cellfun(@(x) or(x, circshift(x, -1)), ContactBeads_index, 'UniformOutput', false);
    ContactBeads_index_aligned = cellfun(@(x) x(1:end-1), ContactBeads_index_aligned, 'UniformOutput', false);

    f_tension_NonContact_magnitude = cellfun(@(x, y) x(~y, :), f_tension_magnitude, ContactBeads_index_aligned, 'UniformOutput', false);
    f_tension_XY = cellfun(@(x) 0.5*(x(1:end-1, :)+x(2:end, :)), XY, 'UniformOutput', false);
    f_tension_XY_NonContact = cellfun(@(x, y) x(~y, :), f_tension_XY, ContactBeads_index_aligned, 'UniformOutput', false);

    % % figure('Position', [100, 100, 800, 800], 'Color', 'w');
    % % for foo = 1:100
    % %     XY_current_in_lattice = mod(XY{foo}-C2C_pillar/2, C2C_pillar);
    % %     ContactBeads_index_current = ContactBeads_index{foo};
    % %     f_tension_XY_NonContact_current = mod(f_tension_XY_NonContact{foo}-C2C_pillar/2, C2C_pillar);
    % %     f_tension_NonContact_magnitude_current = f_tension_NonContact_magnitude{foo};
    % % 
    % %     % plot XY_current_in_lattice in the unit cell with the pillars, each dot represents a bead with a radius of R_fib
    % %     % figure('Position', [100, 100, 800, 800], 'Color', 'w');
    % %     for foo = 1:size(XY_current_in_lattice, 1)
    % %         if ContactBeads_index_current(foo)
    % %             viscircles(XY_current_in_lattice(foo, :)/C2C_pillar, R_fib/C2C_pillar, 'Color', 'm', 'LineStyle', '-'); hold on;
    % %         else
    % %             viscircles(XY_current_in_lattice(foo, :)/C2C_pillar, R_fib/C2C_pillar, 'Color', 'b', 'LineStyle', '-'); hold on;
    % %         end
    % %     end
    % %     viscircles(XY_current_in_lattice(1, :)/C2C_pillar, R_fib/C2C_pillar, 'Color', 'r', 'LineStyle', '-'); hold on;
    % %     viscircles([0.5 0.5], R_pillar/C2C_pillar, 'Color', 'r', 'LineStyle', '--');
    % %     axis equal; axis off
    % %     xlim([0 1]); ylim([0 1]);
    % %     hold on
    % %     scatter(f_tension_XY_NonContact_current(:,1)/C2C_pillar, f_tension_XY_NonContact_current(:,2)/C2C_pillar, 100, f_tension_NonContact_magnitude_current, 'filled');
    % %     colormap('jet'); colorbar;
    % % end

    % Plot the normalized fiber tension in the unit cell
    
    % % % This is the normalized tension by the number of beads under the tension
    % % f_tension_NonContact_magnitude = cellfun(@(x) x/(numel(x)+1), f_tension_NonContact_magnitude, 'UniformOutput', false);

    f_tension_NonContact_magnitude_all = cell2mat(f_tension_NonContact_magnitude);
    f_tension_NonContact_magnitude_all_normalized = f_tension_NonContact_magnitude_all / f_bead;
    f_tension_XY_all = cell2mat(f_tension_XY_NonContact);
    f_tension_XY_in_lattice = mod(f_tension_XY_all-C2C_pillar/2, C2C_pillar);
    f_tension_XY_in_lattice_normalized = f_tension_XY_in_lattice ./ [C2C_pillar, C2C_pillar];
    f_tension_NonContact_normalized_magnitude_in_grid_mean = zeros(length(Xedges)-1, length(Yedges)-1);
    for grid_x = 1:length(Xedges)-1
        for grid_y = 1:length(Yedges)-1
            If_fiber_in_grid = f_tension_XY_in_lattice_normalized(:,1) >= Xedges(grid_x) & ...
                f_tension_XY_in_lattice_normalized(:,1) < Xedges(grid_x+1) & ...
                f_tension_XY_in_lattice_normalized(:,2) >= Yedges(grid_y) & ...
                f_tension_XY_in_lattice_normalized(:,2) < Yedges(grid_y+1);
            f_tension_NonContact_normalized_magnitude_in_grid = f_tension_NonContact_magnitude_all_normalized(If_fiber_in_grid);
            f_tension_NonContact_normalized_magnitude_in_grid_mean(grid_x, grid_y) = mean(f_tension_NonContact_normalized_magnitude_in_grid);
        end
    end
    [XX, YY] = meshgrid(movmean(Xedges, 2, "Endpoints","discard"), movmean(Yedges, 2, "Endpoints","discard"));

    figure('Position', [100, 100, 800, 800], 'Color', 'w');
    pcolor(XX, YY, f_tension_NonContact_normalized_magnitude_in_grid_mean', "EdgeColor", "none");
    % Notice the transpose of the matrix!!!
    axis equal; axis off; grid off
    viscircles([0.5 0.5], R_pillar/C2C_pillar, 'Color', 'r', 'LineStyle', '--');
    xlim([0 1]); ylim([0 1]);
    cmocean('thermal'); colorbar; clim([0 0.2]);
    title('Non-contact tension (normalized), L = ' + string(L_0) + ' um');
    saveas(gcf, [figure_mother_save_path, filesep, 'Non-contact tension (normalized), L = ', L_0{1}, ' um - All.png']);

    L_0_group(ii) = str2double(L_0{1});
    f_tension_NonContact_normalized_magnitude_in_grid_mean_group(:,:,ii) = f_tension_NonContact_normalized_magnitude_in_grid_mean;

    clearvars XY f_tension ContactBeads_index
    close all

end
save([data_mother_save_path, filesep, 'Tension_In_Unit_Cell.mat'], ...
        'L_0_group', 'f_tension_NonContact_normalized_magnitude_in_grid_mean_group', ...
        'XX', 'YY');


%% Load the data in the unit cell, and plot the cut, max. and std. of the normalized tension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

data_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Simulation\Simulations_LubricationModel'];
figure_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Simulation\Tension analysis'];

% load([data_mother_save_path, filesep, 'Tension_In_Unit_Cell.mat']);
load([data_mother_save_path, filesep, 'Tension_In_Unit_Cell_Normalized_by_Bead_Number.mat']);

% Plot a cut of the normalized fiber tension in the unit cell for all L_0
figure('Position', [100, 100, 1000, 600], 'Color', 'w');
cut_position = 43:45;
Vic_cmap = cmocean('amp');
% assign color to each L_0
color_group = Vic_cmap(round(255*(linspace(0, 1, length(L_0_group))))+1, :);
% Order the L_0_group and f_tension_NonContact_normalized_magnitude_in_grid_mean_group based on L_0
[L_0_group, order] = sort(L_0_group);
f_tension_NonContact_normalized_magnitude_in_grid_mean_group = f_tension_NonContact_normalized_magnitude_in_grid_mean_group(:,:,order);
max_tension = zeros(1, length(L_0_group)); 
std_tension = zeros(1, length(L_0_group));
for ii = 1:length(L_0_group)
    % figure('Position', [100, 100, 800, 600], 'Color', 'w');
    plot(YY(:, 1), mean(f_tension_NonContact_normalized_magnitude_in_grid_mean_group(cut_position,:,ii), 1), ...
        'LineWidth', 2, 'Color', color_group(ii, :), 'DisplayName', ['L = ', num2str(L_0_group(ii)), ' um']); 
    hold on
    max_tension(ii) = max(mean(f_tension_NonContact_normalized_magnitude_in_grid_mean_group(cut_position,:,ii), 1));
    std_tension(ii) = std(mean(f_tension_NonContact_normalized_magnitude_in_grid_mean_group(cut_position,:,ii), 1));
end
xlabel('Y position in the unit cell', 'FontSize', 20, 'FontName', 'Times New Roman'); 
ylabel('Normalized tension', 'FontSize', 20, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
legend('Location', 'eastoutside');
title(['Cut of normalized tension (at x = ', num2str(mean(cut_position)/size(YY, 1)), ...
    ') in the unit cell'], 'FontSize', 20);
saveas(gcf, [figure_mother_save_path, filesep, 'Cut of normalized tension in the unit cell (normalized by bead number).png']);

% Plot the max. normalized tension vs. fiber length
figure('Position', [100, 100, 800, 600], 'Color', 'w');
plot(L_0_group, max_tension, 'ko', 'MarkerFaceColor', 'k');
xlabel('Fiber length ($\mu$m)', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('Max. normalized tension', 'FontSize', 20, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
saveas(gcf, [figure_mother_save_path, filesep, 'Max. normalized tension vs. Fiber length (normalized by bead number).png']);

% Plot the std. of normalized tension vs. fiber length
figure('Position', [100, 100, 800, 600], 'Color', 'w');
plot(L_0_group, std_tension, 'ro', 'MarkerFaceColor', 'r');
xlabel('Fiber length ($\mu$m)', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('Std. of normalized tension', 'FontSize', 20, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
saveas(gcf, [figure_mother_save_path, filesep, 'Std. of normalized tension vs. Fiber length (normalized by bead number).png']);




%% Load saved tension data, and max. tension vs x-position for each trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

% Geometrical properties of the array
R_pillar = 8.3 * 1e-6; % pillar radius, unit: m
C2C_pillar = 30 * 1e-6; % m (delta_x & delta_y)
delta_t = 0.1; % s

R_fib = 3e-7; % fiber radius, unit: m
u_max = 200e-6;
mu = 6.1e-3;
f_bead = 6*pi*mu*R_fib*u_max;

Xedges = 0:0.02:1; Yedges = 0:0.02:1; % for the histogram

data_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Simulation\Simulations_LubricationModel'];
saved_data = dir(data_mother_save_path);
% only keep the .mat files inluding the '_TensionStudy'
saved_data = saved_data(arrayfun(@(x) contains(x.name, '_TensionStudy.mat'), saved_data)); 

figure_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Simulation\Tension analysis'];

for ii = 1:length(saved_data)

    load(fullfile(saved_data(ii).folder, saved_data(ii).name));

    L_0 = extractBetween(saved_data(ii).name, 'FiberLength=', 'um_TensionStudy.mat'); % unit: um
    bead_num = size(XY{1, 1}, 1);

    % Apply cellfun to shift contact beads index to align with the tension data
    % Exclude the tension calculated for beads that are in contact with the pillars
    ContactBeads_index_aligned = cellfun(@(x) or(x, circshift(x, -1)), ContactBeads_index, 'UniformOutput', false);
    ContactBeads_index_aligned = cellfun(@(x) x(1:end-1), ContactBeads_index_aligned, 'UniformOutput', false);

    If_contact_all = cellfun(@(x) any(x), ContactBeads_index_aligned);

    f_tension_NonContact_magnitude = cellfun(@(x, y) x(~y, :), f_tension_magnitude, ContactBeads_index_aligned, 'UniformOutput', false);
    f_tension_NonContact_magnitude = cellfun(@(x, y) x/sum(~y), f_tension_NonContact_magnitude, ContactBeads_index_aligned, 'UniformOutput', false);
    f_tension_NonContact_magnitude_max = cellfun(@(x) max(x), f_tension_NonContact_magnitude, 'UniformOutput', false);

    f_tension_XY = cellfun(@(x) 0.5*(x(1:end-1, :)+x(2:end, :)), XY, 'UniformOutput', false);
    f_tension_XY_NonContact = cellfun(@(x, y) x(~y, :), f_tension_XY, ContactBeads_index_aligned, 'UniformOutput', false);
    f_tension_NonContact_magnitude_max_XY = cellfun(@(x, y, z) x(y == z, :), ...
        f_tension_XY_NonContact, f_tension_NonContact_magnitude, f_tension_NonContact_magnitude_max, 'UniformOutput', false);

    % remove the empty cell (notice the order of the code)
    If_contact_all = If_contact_all(~cellfun('isempty', f_tension_NonContact_magnitude_max));
    f_tension_NonContact_magnitude_max_XY = f_tension_NonContact_magnitude_max_XY(~cellfun('isempty', f_tension_NonContact_magnitude_max));
    f_tension_NonContact_magnitude_max = f_tension_NonContact_magnitude_max(~cellfun('isempty', f_tension_NonContact_magnitude_max));

    % only keep the first row in each cell
    f_tension_NonContact_magnitude_max_XY = cellfun(@(x) x(1, :), f_tension_NonContact_magnitude_max_XY, 'UniformOutput', false);

    f_tension_NonContact_magnitude_max_all = cell2mat(f_tension_NonContact_magnitude_max);
    f_tension_NonContact_magnitude_max_XY_all = cell2mat(f_tension_NonContact_magnitude_max_XY);

    x_position = f_tension_NonContact_magnitude_max_XY_all(:, 1);
    % determine the trajectories based on the x-position by finding the sudden drops in the x-position
    x_diff = diff(x_position);
    sudden_drops = find(x_diff < -10*C2C_pillar); % threshold for sudden drop
    start_idx = 1;
    for jj = 1:length(sudden_drops)

        end_idx = sudden_drops(jj);

        current_trajectory_x = f_tension_NonContact_magnitude_max_XY_all(start_idx:end_idx, 1);
        current_trajectory_tension_max = f_tension_NonContact_magnitude_max_all(start_idx:end_idx, 1);
        current_trajectory_If_contact = If_contact_all(start_idx:end_idx, 1);

        Plot_Start_x_indice = sum(current_trajectory_x/C2C_pillar < 100); 
        Plot_Stop_x_indice = sum(current_trajectory_x/C2C_pillar < 150); 
        Plot_trajectory_x = current_trajectory_x(Plot_Start_x_indice:Plot_Stop_x_indice);
        Plot_trajectory_tension_max = current_trajectory_tension_max(Plot_Start_x_indice:Plot_Stop_x_indice);
        Plot_trajectory_If_contact = current_trajectory_If_contact(Plot_Start_x_indice:Plot_Stop_x_indice);

        figure('Position', [100, 100, 1500, 300], 'Color', 'w');
        % plot the max. tension vs x-position, color the contact case in red and non-contact case in blue (NOT in unit cell)
        plot(Plot_trajectory_x(Plot_trajectory_If_contact, 1)/C2C_pillar, ...
            Plot_trajectory_tension_max(Plot_trajectory_If_contact)/f_bead, 'ro', 'MarkerFaceColor', 'r'); hold on;
        plot(Plot_trajectory_x(~Plot_trajectory_If_contact, 1)/C2C_pillar, ...
            Plot_trajectory_tension_max(~Plot_trajectory_If_contact)/f_bead, 'bo', 'MarkerFaceColor', 'b');
        xlabel('$x/\lambda$', 'FontSize', 20, 'Interpreter','latex');
        ylabel('$F_{\rm{t, max}}/F_{\rm{bead}}$', 'FontSize', 20, 'Interpreter','latex');
        set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
        xlim([100 150]); ylim([0 1]);
        title(['Max. tension vs. x-position, L = ', L_0{1}, ' um'], 'FontSize', 20);
        legend('Contact case', 'Non-contact case', 'Location', 'best');
        saveas(gcf, [figure_mother_save_path, filesep, 'Max tension vs X (normalized by bead number)\Max. tension vs. x-position, L = ', L_0{1}, ' um_', num2str(jj), '.png']);

        close all
        start_idx = end_idx + 1;
    end

    % % % for foo = 2000:2050
    % % %     XY_current_in_lattice = mod(XY{foo}-C2C_pillar/2, C2C_pillar);
    % % %     ContactBeads_index_current = ContactBeads_index{foo};
    % % %     f_tension_XY_NonContact_current = mod(f_tension_XY_NonContact{foo}-C2C_pillar/2, C2C_pillar);
    % % %     f_tension_NonContact_magnitude_current = f_tension_NonContact_magnitude{foo};
    % % % 
    % % %     % plot XY_current_in_lattice in the unit cell with the pillars, each dot represents a bead with a radius of R_fib
    % % %     figure('Position', [100, 100, 800, 800], 'Color', 'w');
    % % %     for fooOO = 1:size(XY_current_in_lattice, 1)
    % % %         if ContactBeads_index_current(fooOO)
    % % %             viscircles(XY_current_in_lattice(fooOO, :)/C2C_pillar, R_fib/C2C_pillar, 'Color', 'm', 'LineStyle', '-'); hold on;
    % % %         else
    % % %             viscircles(XY_current_in_lattice(fooOO, :)/C2C_pillar, R_fib/C2C_pillar, 'Color', 'b', 'LineStyle', '-'); hold on;
    % % %         end
    % % %     end
    % % %     viscircles(XY_current_in_lattice(1, :)/C2C_pillar, R_fib/C2C_pillar, 'Color', 'r', 'LineStyle', '-'); hold on;
    % % %     viscircles([0.5 0.5], R_pillar/C2C_pillar, 'Color', 'r', 'LineStyle', '--');
    % % %     axis equal; % axis off 
    % % %     % xlim([0 1]); ylim([0 1]);
    % % %     hold on
    % % %     scatter(f_tension_XY_NonContact_current(:,1)/C2C_pillar, f_tension_XY_NonContact_current(:,2)/C2C_pillar, ...
    % % %         60, f_tension_NonContact_magnitude_current/f_bead, 'filled');
    % % %     colormap('jet'); colorbar;
    % % %     set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
    % % % 
    % % %     close all;
    % % % end

    clearvars XY f_tension ContactBeads_index
end



%% Load saved tension data, plot ensemble average results for each length, and best metric parameter vs. fiber length (Fig.5 in the paper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

% Geometrical properties of the array
R_pillar = 8.3 * 1e-6; % pillar radius, unit: m
C2C_pillar = 30 * 1e-6; % m (delta_x & delta_y)
delta_t = 0.1; % s

R_fib = 3e-7; % fiber radius, unit: m
u_max = 200e-6;
mu = 6.1e-3;
f_bead = 6*pi*mu*R_fib*u_max;

Xedges = 0:0.02:1; Yedges = 0:0.02:1; % for the histogram

data_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Simulation\Simulations_LubricationModel'];
saved_data = dir(data_mother_save_path);
% only keep the .mat files inluding the '_TensionStudy'
saved_data = saved_data(arrayfun(@(x) contains(x.name, '_TensionStudy.mat'), saved_data)); 

figure_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Simulation\Tension analysis'];

L_0_group = zeros(1, length(saved_data));
max_ensemble_tension = zeros(1, length(saved_data));
max_ensemble_tension_std = zeros(1, length(saved_data));
for ii = 1:length(saved_data)

    load(fullfile(saved_data(ii).folder, saved_data(ii).name));

    L_0 = extractBetween(saved_data(ii).name, 'FiberLength=', 'um_TensionStudy.mat'); % unit: um
    bead_num = size(XY{1, 1}, 1);

    f_tension_XY = cellfun(@(x) 0.5*(x(1:end-1, :)+x(2:end, :)), XY, 'UniformOutput', false);

    %%%%%%%%%%%%%%%%%%%%% To select the contact cases %%%%%%%%%%%%%%%%%%%%%
    If_contact_all = cellfun(@(x) any(x), ContactBeads_index);
    % Choose only the contact cases
    f_tension_magnitude_Contact = f_tension_magnitude(If_contact_all);
    f_tension_XY_Contact = f_tension_XY(If_contact_all);
    ContactBeads_index_Contact = ContactBeads_index(If_contact_all);
    %%%%%%%%%%%%%%%%%%%%% To select the contact cases %%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%% To check the order of the beads in the fiber %%%%%%%%%%%%%%%%%%%%%
    theta = 35; % angle in degrees (counterclockwise)
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; % rotation matrix

    f_tension_XY_Contact_rotated = cellfun(@(x) (R * x')', f_tension_XY_Contact, 'UniformOutput', false);
    f_tension_X_Contact_rotated = cellfun(@(x) x(:, 1), f_tension_XY_Contact_rotated, 'UniformOutput', false);

    beads_order_indicator = cellfun(@(x) x(end)-x(1), f_tension_X_Contact_rotated, 'UniformOutput', false);
    If_inOrder = cellfun(@(x) x>=0, beads_order_indicator, 'UniformOutput', false);

    % Plot the rotated fiber (color indicates the order of the beads in the fiber)
    % % % for foo = 1:30:1000
    % % %     XY_current = f_tension_XY_Contact_rotated{foo};
    % % %     If_inOrder_current = If_inOrder{foo};
    % % % 
    % % %     % plot XY_current in the unit cell with the pillars, each dot represents a bead with a radius of R_fib
    % % %     figure('Position', [100, 100, 800, 600], 'Color', 'w');
    % % %     if If_inOrder_current
    % % %         plot(XY_current(:, 1), XY_current(:, 2), 'b', 'LineWidth', 5); hold on;
    % % %     else
    % % %         plot(XY_current(:, 1), XY_current(:, 2), 'r', 'LineWidth', 5); hold on;
    % % %     end            
    % % %     viscircles(XY_current(1, :), R_fib*2, 'Color', 'm', 'LineStyle', '-'); hold on;
    % % %     axis equal; axis off 
    % % %     set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
    % % % end
    %%%%%%%%%%%%%%%%%%%%% To check the order of the beads in the fiber %%%%%%%%%%%%%%%%%%%%%

    If_inOrder_all = cell2mat(If_inOrder);
    
    f_tension_magnitude_Contact_all = reshape(cell2mat(f_tension_magnitude_Contact), [], numel(f_tension_magnitude_Contact));
    ContactBeads_index_Contact_all = reshape(cell2mat(ContactBeads_index_Contact), [], numel(ContactBeads_index_Contact));
    ContactBeads_index_Contact_all(end, :) = []; % remove the last row (to align with the f_tension_magnitude_Contact_all)

    % Flip the order of f_tension_magnitude_Contact_all based on If_inOrder_all
    f_tension_magnitude_Contact_all(:, ~If_inOrder_all) = flipud(f_tension_magnitude_Contact_all(:, ~If_inOrder_all));
    ContactBeads_index_Contact_all(:, ~If_inOrder_all) = flipud(ContactBeads_index_Contact_all(:, ~If_inOrder_all));

    % set the f_tension_magnitude to NaN for the beads that are in contact with the pillars
    f_tension_magnitude_Contact_all(ContactBeads_index_Contact_all == 1) = NaN;

    f_tension_magnitude_Contact_average = mean(f_tension_magnitude_Contact_all, 2, 'omitnan');
    f_tension_magnitude_Contact_std = std(f_tension_magnitude_Contact_all, 0, 2, 'omitnan');

    figure('Position', [100, 100, 800, 600], 'Color', 'w');
    plot((1:numel(f_tension_magnitude_Contact_average))'/bead_num, f_tension_magnitude_Contact_average/f_bead, 'k--', 'LineWidth', 3); hold on;
    for jj = 1:round(size(f_tension_magnitude_Contact_all, 2)/100):size(f_tension_magnitude_Contact_all, 2)
        plot((1:numel(f_tension_magnitude_Contact_average))'/bead_num, f_tension_magnitude_Contact_all(:, jj)/f_bead, 'Color', [1 0 1 0.1]); hold on;
    end
    plot((1:numel(f_tension_magnitude_Contact_average))'/bead_num, f_tension_magnitude_Contact_average/f_bead, 'k--', 'LineWidth', 3); hold on;
    xlim([0 1]); ylim([0 8]);
    legend('Average', 'Location', 'northeast', 'FontSize', 40, 'FontName', 'Times New Roman');
    xlabel('$s/L$', 'FontSize', 40, 'Interpreter', 'latex');
    ylabel('$\tilde{F}_{\rm{tension}}$', 'FontSize', 40, 'Interpreter', 'latex');
    % title(['Tension along the fiber, L = ', L_0{1}, ' um'], 'FontSize', 20);
    text(0.25, 7.1, ['L = ', L_0{1}, ' ', char(181) 'm'], 'FontSize', 40, 'FontName', 'Times New Roman', 'HorizontalAlignment', 'center');
    set(gca, 'FontSize', 40, 'FontName', 'Times New Roman');
    grid on
    set(gcf,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',[figure_mother_save_path, filesep, ...
        'Normalized tension along the fiber, only contact cases and beads reordered downstream, L = ', L_0{1}, ' um.eps']);
    % saveas(gcf, [figure_mother_save_path, filesep, 'Normalized tension along the fiber, ' ...
    %     'only contact cases and beads reordered downstream, L = ', L_0{1}, ' um.png']);
  
    f_tension_X_Contact_rotated_all = reshape(cell2mat(f_tension_X_Contact_rotated), [], numel(f_tension_X_Contact_rotated));
    x_position = mean(f_tension_X_Contact_rotated_all);
    % determine the trajectories based on the x-position by finding the sudden drops in the x-position
    x_diff = diff(x_position);
    sudden_drops = find(x_diff < -10*C2C_pillar); % threshold for sudden drop
    start_idx = 1;
    max_ensemble_tension_per_case = zeros(numel(sudden_drops), 1);
    for jj = 1:length(sudden_drops)

        end_idx = sudden_drops(jj);

        current_f_tension_magnitude = f_tension_magnitude_Contact_all(:, start_idx:end_idx);
        current_f_tension_magnitude_average = mean(current_f_tension_magnitude, 2, 'omitnan');

        max_ensemble_tension_per_case(jj) = max(current_f_tension_magnitude_average)/f_bead;
        
        start_idx = end_idx + 1;
    end

    L_0_group(ii) = str2double(L_0{1});
    max_ensemble_tension(ii) = mean(max_ensemble_tension_per_case);
    max_ensemble_tension_std(ii) = std(max_ensemble_tension_per_case);

    close all
    clearvars XY f_tension ContactBeads_index
end

% Plot the max. normalized tension vs. fiber length
figure('Position', [100, 100, 800, 600], 'Color', 'w');
yyaxis left
ax = gca; 
ax.YColor = [231, 138, 195]/256;
errorbar(L_0_group*1e-6/C2C_pillar, max_ensemble_tension, max_ensemble_tension_std, 'o', ...
     'MarkerFaceColor', [231, 138, 195]/256, 'MarkerEdgeColor', 'k', ...
     'MarkerSize', 15, 'LineWidth', 1.5, 'Color', 'k');
xlabel('$L/\lambda$', 'FontSize', 28, 'Interpreter', 'latex');
ylabel('$\overline{\langle \tilde{F}_{\rm{tension,\,contact}} \rangle {_{\rm{max}}}}$', 'FontSize', 28, 'Interpreter', 'latex');
xlim([0 3.2]); ylim([0 5]); % grid on
set(gca, 'FontSize', 28, 'FontName', 'Times New Roman');

% See the following code: VicFig_Get_fig_data.m

% % saveas(gcf, [figure_mother_save_path, filesep, 'Max. ensemble averaged normalized tension vs. Fiber length.png']);




%% Load saved tension data, plot chronophotography of two cases in regimes I and regime II, and their tension along the fiber (Fig.3 in the paper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

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

X_pillar = 0:C2C_pillar:1e-2; Y_pillar = 0:-C2C_pillar:-1e-2;
[X_pillar, Y_pillar] = meshgrid(X_pillar, Y_pillar);

theta = -35; % angle in degrees (clockwise)
Rotate_matrix = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; % rotation matrix

XY_pillar_rotated = [X_pillar(:), Y_pillar(:)]*Rotate_matrix;
% % figure('Position', [100, 100, 800, 800], 'Color', 'w');
% % hold on;
% % viscircles(XY_pillar_rotated, R_pillar*ones(size(XY_pillar_rotated, 1), 1), 'Color', 'r');
% % axis equal; axis off

figure_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Simulation\Tension analysis'];

parent_path = 'E:\Experimental Data (RAW)\Simulation\data';
sub1_path = dir(parent_path);

% % % % Case: n30, no.17 (Short fiber)
% % % sub1Path_i = 11; sub2Path_i = 14;
% Case: n90, no.1 (Long fiber)
sub1Path_i = 17; sub2Path_i = 6;

current_beadsNum = sub1_path(sub1Path_i).name;
beads_num = str2double(current_beadsNum(2:end));
disp(['Current beads number: ', num2str(beads_num)]);

L_0 = 0.6 * beads_num * 1e-6;
disp(['Current fiber length: ', num2str(L_0), ' um']);

sub2_path = dir(fullfile(parent_path, current_beadsNum));

current_sim = sub2_path(sub2Path_i).name;
sim_no = str2double(current_sim(11:end));
disp(['Current simulation number: ', num2str(sim_no)]);

fileinfo = dir(fullfile(parent_path, current_beadsNum, current_sim, 'output_data\*.vtk'));

XY = zeros(length(fileinfo), beads_num, 2);
ContactBeads_index_current = zeros(length(fileinfo), beads_num);
XY_tension = zeros(length(fileinfo), beads_num-1, 2);
XY_tension_rotated = zeros(length(fileinfo), beads_num-1, 2);
f_tension_magnitude = zeros(length(fileinfo), beads_num-1);
for ii = 1:length(fileinfo)
    snapshot = readVTK(fullfile(fileinfo(ii).folder, fileinfo(ii).name));
    XY_current = snapshot.points(:, 1:2);

    XY_current_in_lattice = mod(XY_current-C2C_pillar/2, C2C_pillar);
    dist2origin = sqrt((XY_current_in_lattice(:,1)-C2C_pillar/2).^2 + (XY_current_in_lattice(:,2)-C2C_pillar/2).^2);
    ContactBeads_index_current(ii, :) = dist2origin < R_pillar + 0.3 * 1e-6;

    XY(ii, :, :) = XY_current;
    XY_tension(ii, :, :) = 0.5*(XY_current(1:end-1, :)+XY_current(2:end, :));
    XY_tension_rotated(ii, :, :) = 0.5*(XY_current(1:end-1, :)+XY_current(2:end, :))*Rotate_matrix;

    vec_im = diff(XY_current);
    l_im = vecnorm(vec_im')';
    vec_t = vec_im./l_im;
    f_tension_magnitude(ii, :) = ks*(l_im - 2*R_fib); % Notice the sign
end

% % % %%%%%% Plot the tension along the fiber (Short fiber) %%%%%%
% % % figure('Position', [100, 100, 800, 600], 'Color', 'w');
% % % unused_color_ratio = 0.4;
% % % 
% % % start_snapshot = 49; end_snapshot = 68; 
% % % picked_snaps = start_snapshot:2:end_snapshot;
% % % logical(sum(ContactBeads_index_current(picked_snaps, :), 2));
% % % 
% % % picked_snaps_noncontact = picked_snaps(~logical(sum(ContactBeads_index_current(picked_snaps, :), 2)));
% % % total_snapshots = numel(picked_snaps_noncontact);
% % % color_group_noncontact = cmocean('algae', round(total_snapshots/(1-unused_color_ratio)));
% % % f_tension_magnitude_plotted = zeros(numel(picked_snaps), beads_num-1);
% % % n = 1; no_color = 1; 
% % % for ii = picked_snaps_noncontact
% % %     f_tension_magnitude_to_plot = f_tension_magnitude(ii, :)/f_bead;
% % %     if XY(ii, 1, 1) > XY(ii, end, 1)
% % %         f_tension_magnitude_to_plot = flip(f_tension_magnitude_to_plot);
% % %     end
% % %     plot((1:beads_num-1)/beads_num, f_tension_magnitude_to_plot, 'Color', ...
% % %         color_group_noncontact(no_color+round(total_snapshots*unused_color_ratio/2), :), 'LineWidth', 1.5); 
% % %     hold on;
% % %     f_tension_magnitude_plotted(n, :) = f_tension_magnitude_to_plot;
% % %     n = n + 1; no_color = no_color + 1;
% % % end
% % % 
% % % picked_snaps_contact = picked_snaps(logical(sum(ContactBeads_index_current(picked_snaps, :), 2)));
% % % total_snapshots = numel(picked_snaps_contact);
% % % color_group_contact = cmocean('solar', round(total_snapshots/(1-unused_color_ratio)));
% % % no_color = 1;
% % % for ii = picked_snaps_contact
% % %     f_tension_magnitude_to_plot = f_tension_magnitude(ii, :)/f_bead;
% % %     f_tension_magnitude_to_plot(logical(ContactBeads_index_current(ii, 1:end-1))) = NaN;
% % %     if XY(ii, 1, 1) > XY(ii, end, 1)
% % %         f_tension_magnitude_to_plot = flip(f_tension_magnitude_to_plot);
% % %     end
% % %     plot((1:beads_num-1)/beads_num, f_tension_magnitude_to_plot, 'Color', ...
% % %         color_group_contact(no_color+round(total_snapshots*unused_color_ratio/2), :), 'LineWidth', 1.5); 
% % %     hold on;
% % %     f_tension_magnitude_plotted(n, :) = f_tension_magnitude_to_plot;
% % %     n = n + 1; no_color = no_color + 1;
% % % end
% % % plot((1:beads_num-1)/beads_num, mean(f_tension_magnitude_plotted, 1, 'omitnan'), 'm-', 'LineWidth', 3);
% % % xlabel('$s/L$', 'FontSize', 28, 'Interpreter', 'latex');
% % % ylabel('$\tilde{F}_{\rm{tension}}$', 'FontSize', 28, 'Interpreter', 'latex');
% % % xlim([0 1]); ylim([0 8]); grid on
% % % set(gca, 'FontSize', 28, 'FontName', 'Times New Roman');
% % % set(gcf,'renderer','Painters');
% % % print('-depsc2','-tiff','-r100','-vector',[figure_mother_save_path, filesep, 'Short_fiber_tension_profile_new2.eps']);
% % % 
% % % % Plot the chronophotography (Short fiber)
% % % figure('Position', [100, 100, 800, 600], 'Color', 'w');
% % % no_color = 1;
% % % for ii = picked_snaps_noncontact
% % %     plot(XY_tension_rotated(ii, :, 1), XY_tension_rotated(ii, :, 2), 'Color', ...
% % %         color_group_noncontact(no_color+round(total_snapshots*unused_color_ratio/2), :), 'LineWidth', 2);
% % %     hold on;
% % %     no_color = no_color + 1;
% % % end
% % % no_color = 1;
% % % for ii = picked_snaps_contact
% % %     plot(XY_tension_rotated(ii, :, 1), XY_tension_rotated(ii, :, 2), 'Color', ...
% % %         color_group_contact(no_color+round(total_snapshots*unused_color_ratio/2), :), 'LineWidth', 2);
% % %     hold on;
% % %     no_color = no_color + 1;
% % % end
% % % XY_pillar_rotated(XY_pillar_rotated(:, 2) < -5e-5, :) = [];
% % % XY_pillar_rotated(XY_pillar_rotated(:, 2) > 8e-5, :) = [];
% % % XY_pillar_rotated(XY_pillar_rotated(:, 1) < 4e-4, :) = [];
% % % XY_pillar_rotated(XY_pillar_rotated(:, 1) > 7e-4, :) = [];
% % % viscircles(XY_pillar_rotated, R_pillar*ones(size(XY_pillar_rotated, 1), 1), 'Color', 'k', 'LineWidth', 0.2);
% % % axis equal; axis off
% % % xlim([4.5 6.5]*1e-4); ylim([0.5 4.5]*1e-5);
% % % set(gcf,'renderer','Painters');
% % % print('-depsc2','-tiff','-r100','-vector',[figure_mother_save_path, filesep, 'Short_fiber_chronophotography_new2.eps']);


%%%%%% Plot the tension along the fiber (Long fiber) %%%%%%
figure('Position', [100, 100, 800, 600], 'Color', 'w');
unused_color_ratio = 0.4;

start_snapshot = 50; end_snapshot = 69; 
picked_snaps = start_snapshot:2:end_snapshot;
logical(sum(ContactBeads_index_current(picked_snaps, :), 2));

picked_snaps_noncontact = picked_snaps(~logical(sum(ContactBeads_index_current(picked_snaps, :), 2)));
total_snapshots = numel(picked_snaps_noncontact);
color_group_noncontact = cmocean('algae', round(total_snapshots/(1-unused_color_ratio)));
f_tension_magnitude_plotted = zeros(numel(picked_snaps), beads_num-1);
n = 1; no_color = 1; 
for ii = picked_snaps_noncontact
    f_tension_magnitude_to_plot = f_tension_magnitude(ii, :)/f_bead;
    if XY(ii, 1, 1) > XY(ii, end, 1) 
        f_tension_magnitude_to_plot = flip(f_tension_magnitude_to_plot);
    end
    plot((1:beads_num-1)/beads_num, f_tension_magnitude_to_plot, 'Color', ...
        color_group_noncontact(no_color+round(total_snapshots*unused_color_ratio/2), :), 'LineWidth', 1.5); 
    hold on;
    f_tension_magnitude_plotted(n, :) = f_tension_magnitude_to_plot;
    n = n + 1; no_color = no_color + 1;
end

picked_snaps_contact = picked_snaps(logical(sum(ContactBeads_index_current(picked_snaps, :), 2)));
total_snapshots = numel(picked_snaps_contact);
color_group_contact = cmocean('solar', round(total_snapshots/(1-unused_color_ratio)));
no_color = 1;
for ii = picked_snaps_contact
    f_tension_magnitude_to_plot = f_tension_magnitude(ii, :)/f_bead;
    f_tension_magnitude_to_plot(logical(ContactBeads_index_current(ii, 1:end-1))) = NaN;
    if XY(ii, 1, 1) > XY(ii, end, 1) 
        f_tension_magnitude_to_plot = flip(f_tension_magnitude_to_plot);
    end
    plot((1:beads_num-1)/beads_num, f_tension_magnitude_to_plot, 'Color', ...
        color_group_contact(no_color+round(total_snapshots*unused_color_ratio/2), :), 'LineWidth', 1.5); 
    hold on;
    f_tension_magnitude_plotted(n, :) = f_tension_magnitude_to_plot;
    n = n + 1; no_color = no_color + 1;
end
plot((1:beads_num-1)/beads_num, mean(f_tension_magnitude_plotted, 1, 'omitnan'), 'm-', 'LineWidth', 3);
xlabel('$s/L$', 'FontSize', 28, 'Interpreter', 'latex');
ylabel('$\tilde{F}_{\rm{tension}}$', 'FontSize', 28, 'Interpreter', 'latex');
xlim([0 1]); ylim([0 8]); grid on
set(gca, 'FontSize', 28, 'FontName', 'Times New Roman');
set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',[figure_mother_save_path, filesep, 'Long_fiber_tension_profile_new2.eps']);

% Plot the chronophotography (Long fiber)
figure('Position', [100, 100, 800, 600], 'Color', 'w');
no_color = 1;
for ii = picked_snaps_noncontact
    plot(XY_tension_rotated(ii, :, 1), XY_tension_rotated(ii, :, 2), 'Color', ...
        color_group_noncontact(no_color+round(total_snapshots*unused_color_ratio/2), :), 'LineWidth', 2);
    hold on;
    no_color = no_color + 1;
end
no_color = 1;
for ii = picked_snaps_contact
    plot(XY_tension_rotated(ii, :, 1), XY_tension_rotated(ii, :, 2), 'Color', ...
        color_group_contact(no_color+round(total_snapshots*unused_color_ratio/2), :), 'LineWidth', 2);
    hold on;
    no_color = no_color + 1;
end
XY_pillar_rotated(XY_pillar_rotated(:, 2) < 0.2e-4, :) = [];
XY_pillar_rotated(XY_pillar_rotated(:, 2) > 1.3e-4, :) = [];
XY_pillar_rotated(XY_pillar_rotated(:, 1) < 2.5e-4, :) = [];
XY_pillar_rotated(XY_pillar_rotated(:, 1) > 5.5e-4, :) = [];
viscircles(XY_pillar_rotated, R_pillar*ones(size(XY_pillar_rotated, 1), 1), 'Color', 'k', 'LineWidth', 0.2);
axis equal; axis off
xlim([2.7 4.7]*1e-4); ylim([5.5 9.5]*1e-5);
set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',[figure_mother_save_path, filesep, 'Long_fiber_chronophotography_new2.eps']);


% % % Plot the colorbrar
% % % figure('color', 'w'); set(gcf, 'Position', [100 100 1800 400]);
% % % cmocean('tempo', round(total_snapshots/(1-unused_color_ratio)));
% % % clim([1 round(total_snapshots/(1-unused_color_ratio))+1])
% % % time_start = 1+round(total_snapshots*unused_color_ratio/2);
% % % c = colorbar;
% % % c.Label.String = 'Time: $s$';
% % % c.Label.Interpreter = 'LaTeX';
% % % c.TickLabelInterpreter = 'LaTeX';
% % % c.FontSize = 42;
% % % c.Location = 'north';
% % % if mod(time_start, 2) == 0
% % %     c.Ticks = 0:2:round(total_snapshots/(1-unused_color_ratio));
% % %     c.TickLabels = arrayfun(@(x) sprintf('%.1f', x*time_step), ...
% % %     -time_start:2:round(total_snapshots/(1-unused_color_ratio))-time_start, 'UniformOutput', false);
% % % else
% % %     c.Ticks = 1:2:round(total_snapshots/(1-unused_color_ratio));
% % %     c.TickLabels = arrayfun(@(x) sprintf('%.1f', x*time_step), ...
% % %     -time_start+1:2:round(total_snapshots/(1-unused_color_ratio))-time_start+1, 'UniformOutput', false);
% % % end
% % % axis off
% % % set(gcf,'renderer','Painters');
% % % print('-depsc2','-tiff','-r100','-vector',[figure_mother_save_path, filesep, 'Fiber_chronophotography_colorbar.eps']);



%% Load saved tension data, plot fiber length over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

R_fib = 3e-7; % fiber radius, unit: m

parent_path = 'E:\Experimental Data (RAW)\Simulation\data';
sub1_path = dir(parent_path);

% Case: n30, no.17 (Short fiber)
sub1Path_i = 11; sub2Path_i = 14;
% Case: n90, no.1 (Long fiber)
sub1Path_i = 17; sub2Path_i = 6;

current_beadsNum = sub1_path(sub1Path_i).name;
beads_num = str2double(current_beadsNum(2:end));
disp(['Current beads number: ', num2str(beads_num)]);

L_0 = 0.6 * beads_num * 1e-6;
disp(['Current fiber length: ', num2str(L_0), ' um']);

sub2_path = dir(fullfile(parent_path, current_beadsNum));

current_sim = sub2_path(sub2Path_i).name;
sim_no = str2double(current_sim(11:end));
disp(['Current simulation number: ', num2str(sim_no)]);

fileinfo = dir(fullfile(parent_path, current_beadsNum, current_sim, 'output_data\*.vtk'));

L_real = zeros(length(fileinfo), 1);
for ii = 1:length(fileinfo)
    snapshot = readVTK(fullfile(fileinfo(ii).folder, fileinfo(ii).name));
    XY_current = snapshot.points(:, 1:2);

    % % figure('Position', [100, 100, 800, 800], 'Color', 'w');
    % % plot(XY_current(:,1), XY_current(:,2), '.', 'LineWidth', 2);
    % % axis equal

    L_real(ii) = sum(sqrt(sum(diff(XY_current).^2, 2))) + 2*R_fib;
end

% Plot the tension along the fiber (Short fiber)
figure('Position', [100, 100, 800, 600], 'Color', 'w');
plot(L_real / L_0, 'k-', 'LineWidth', 2);
xlabel('Snapshot Index', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('Normalized Fiber Length', 'FontSize', 20, 'FontName', 'Times New Roman');
title('Normalized Fiber Length over Time', 'FontSize', 20, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
grid on;





%% Load saved tension data, plot chronophotography of two cases in regimes I and regime II, and average tension over time (new Fig.3 in the paper)
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

X_pillar = 0:C2C_pillar:1e-2; Y_pillar = 0:-C2C_pillar:-1e-2;
[X_pillar, Y_pillar] = meshgrid(X_pillar, Y_pillar);

theta = -35; % angle in degrees (clockwise)
Rotate_matrix = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; % rotation matrix

XY_pillar_rotated = [X_pillar(:), Y_pillar(:)]*Rotate_matrix;
% % figure('Position', [100, 100, 800, 800], 'Color', 'w');
% % hold on;
% % viscircles(XY_pillar_rotated, R_pillar*ones(size(XY_pillar_rotated, 1), 1), 'Color', 'r');
% % axis equal; axis off

figure_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Simulation\Tension analysis'];

parent_path = 'E:\Experimental Data (RAW)\Simulation\data';
sub1_path = dir(parent_path);

for case_no = 1:2

    if case_no == 1
        % Case: n30, no.17 (Short fiber)
        sub1Path_i = 11; sub2Path_i = 14;
    elseif case_no == 2
        % Case: n90, no.1 (Long fiber)
        sub1Path_i = 17; sub2Path_i = 6;
    end

    current_beadsNum = sub1_path(sub1Path_i).name;
    beads_num = str2double(current_beadsNum(2:end));
    disp(['Current beads number: ', num2str(beads_num)]);

    L_0 = 0.6 * beads_num * 1e-6;
    disp(['Current fiber length: ', num2str(L_0), ' um']);

    sub2_path = dir(fullfile(parent_path, current_beadsNum));

    current_sim = sub2_path(sub2Path_i).name;
    sim_no = str2double(current_sim(11:end));
    disp(['Current simulation number: ', num2str(sim_no)]);

    fileinfo = dir(fullfile(parent_path, current_beadsNum, current_sim, 'output_data\*.vtk'));

    XY = zeros(length(fileinfo), beads_num, 2);
    ContactBeads_index_current = zeros(length(fileinfo), beads_num);
    XY_tension = zeros(length(fileinfo), beads_num-1, 2);
    XY_tension_rotated = zeros(length(fileinfo), beads_num-1, 2);
    f_tension_magnitude = zeros(length(fileinfo), beads_num-1);
    for ii = 1:length(fileinfo)
        snapshot = readVTK(fullfile(fileinfo(ii).folder, fileinfo(ii).name));
        XY_current = snapshot.points(:, 1:2);

        XY_current_in_lattice = mod(XY_current-C2C_pillar/2, C2C_pillar);
        dist2origin = sqrt((XY_current_in_lattice(:,1)-C2C_pillar/2).^2 + (XY_current_in_lattice(:,2)-C2C_pillar/2).^2);
        ContactBeads_index_current(ii, :) = dist2origin < R_pillar + 0.3 * 1e-6;

        XY(ii, :, :) = XY_current;
        XY_tension(ii, :, :) = 0.5*(XY_current(1:end-1, :)+XY_current(2:end, :));
        XY_tension_rotated(ii, :, :) = 0.5*(XY_current(1:end-1, :)+XY_current(2:end, :))*Rotate_matrix;

        vec_im = diff(XY_current);
        l_im = vecnorm(vec_im')';
        vec_t = vec_im./l_im;
        f_tension_magnitude(ii, :) = ks*(l_im - 2*R_fib); % Notice the sign
    end

    markers = {'v', '^'};

    %%%%%% Plot the average tension vs. time (Long and short fiber) %%%%%%
    if case_no ==  1
        figure('Position', [100, 100, 800, 600], 'Color', 'w');
        % ONLY for legend
        plot(NaN, NaN, markers{1}, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', ...
            'k', 'MarkerSize', 15, 'LineWidth', 1.5, 'Color', 'k');
        hold on
        plot(NaN, NaN, markers{2}, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', ...
            'k', 'MarkerSize', 15, 'LineWidth', 1.5, 'Color', 'k');
    end
    unused_color_ratio = 0.4;

    start_snapshot = 55; end_snapshot = 64;
    picked_snaps = start_snapshot:end_snapshot;
    logical(sum(ContactBeads_index_current(picked_snaps, :), 2));

    picked_snaps_noncontact = picked_snaps(~logical(sum(ContactBeads_index_current(picked_snaps, :), 2)));
    total_snapshots = numel(picked_snaps_noncontact);
    color_group_noncontact = cmocean('algae', round(total_snapshots/(1-unused_color_ratio)));
    f_tension_magnitude_plotted = zeros(numel(picked_snaps), beads_num-1);
    n = 1; no_color = 1;
    for ii = picked_snaps_noncontact
        f_tension_magnitude_to_plot = f_tension_magnitude(ii, :)/f_bead;
        f_tension_magnitude_to_plot_mean = mean(f_tension_magnitude_to_plot, 'omitnan');
        plot((ii-picked_snaps(1))*dt, f_tension_magnitude_to_plot_mean, markers{case_no}, ...
            'MarkerFaceColor', color_group_noncontact(no_color+round(total_snapshots*unused_color_ratio/2), :), ...
            'MarkerEdgeColor', 'k', 'MarkerSize', 15, 'LineWidth', 1.5, 'Color', 'k');
        hold on;
        f_tension_magnitude_plotted(n, :) = f_tension_magnitude_to_plot;
        n = n + 1; no_color = no_color + 1;
    end

    picked_snaps_contact = picked_snaps(logical(sum(ContactBeads_index_current(picked_snaps, :), 2)));
    total_snapshots = numel(picked_snaps_contact);
    color_group_contact = cmocean('solar', round(total_snapshots/(1-unused_color_ratio)));
    no_color = 1;
    for ii = picked_snaps_contact
        f_tension_magnitude_to_plot = f_tension_magnitude(ii, :)/f_bead;
        f_tension_magnitude_to_plot(logical(ContactBeads_index_current(ii, 1:end-1))) = NaN;
        f_tension_magnitude_to_plot_mean = mean(f_tension_magnitude_to_plot, 'omitnan');
        plot((ii-picked_snaps(1))*dt, f_tension_magnitude_to_plot_mean, markers{case_no}, ...
            'MarkerFaceColor', color_group_contact(no_color+round(total_snapshots*unused_color_ratio/2), :), ...
            'MarkerEdgeColor', 'k', 'MarkerSize', 15, 'LineWidth', 1.5, 'Color', 'k');
        hold on;
        f_tension_magnitude_plotted(n, :) = f_tension_magnitude_to_plot;
        n = n + 1; no_color = no_color + 1;
    end
    % plot((1:beads_num-1)/beads_num, mean(f_tension_magnitude_plotted, 1, 'omitnan'), 'm-', 'LineWidth', 3);
    xlabel('$\rm{Time\,(s)}$', 'FontSize', 28, 'Interpreter', 'latex');
    ylabel('$\tilde{F}_{\rm{tension,\,mean}}$', 'FontSize', 28, 'Interpreter', 'latex');
    xlim([0 1]); ylim([0 4]); grid on
    set(gca, 'FontSize', 28, 'FontName', 'Times New Roman');
    set(gcf,'renderer','Painters');

    legend({'Short fiber', 'Long fiber'}, 'Location', 'northeast', 'FontSize', 28, 'FontName', 'Times New Roman');
    % print('-depsc2','-tiff','-r100','-vector',[figure_mother_save_path, filesep, 'Long_and_short_fiber_tension_vs_time.eps']);


    % if case_no == 1
    %     % Plot the chronophotography (Short fiber)
    %     figure('Position', [100, 100, 800, 600], 'Color', 'w');
    %     no_color = 1;
    %     for ii = picked_snaps_noncontact
    %         plot(XY_tension_rotated(ii, :, 1), XY_tension_rotated(ii, :, 2), 'Color', ...
    %             color_group_noncontact(no_color+round(total_snapshots*unused_color_ratio/2), :), 'LineWidth', 2);
    %         hold on;
    %         no_color = no_color + 1;
    %     end
    %     no_color = 1;
    %     for ii = picked_snaps_contact
    %         plot(XY_tension_rotated(ii, :, 1), XY_tension_rotated(ii, :, 2), 'Color', ...
    %             color_group_contact(no_color+round(total_snapshots*unused_color_ratio/2), :), 'LineWidth', 2);
    %         hold on;
    %         no_color = no_color + 1;
    %     end
    %     XY_pillar_rotated(XY_pillar_rotated(:, 2) < -5e-5, :) = [];
    %     XY_pillar_rotated(XY_pillar_rotated(:, 2) > 8e-5, :) = [];
    %     XY_pillar_rotated(XY_pillar_rotated(:, 1) < 4e-4, :) = [];
    %     XY_pillar_rotated(XY_pillar_rotated(:, 1) > 7e-4, :) = [];
    %     viscircles(XY_pillar_rotated, R_pillar*ones(size(XY_pillar_rotated, 1), 1), 'Color', 'k', 'LineWidth', 0.2);
    %     axis equal; axis off
    %     xlim([4.8 6.2]*1e-4); ylim([0 5]*1e-5);
    %     set(gcf,'renderer','Painters');
    %     print('-depsc2','-tiff','-r100','-vector',[figure_mother_save_path, filesep, 'Short_fiber_chronophotography_new.eps']);
    % 
    % elseif case_no == 2
    %     % Plot the chronophotography (Long fiber)
    %     figure('Position', [100, 100, 800, 600], 'Color', 'w');
    %     no_color = 1;
    %     for ii = picked_snaps_noncontact
    %         plot(XY_tension_rotated(ii, :, 1), XY_tension_rotated(ii, :, 2), 'Color', ...
    %             color_group_noncontact(no_color+round(total_snapshots*unused_color_ratio/2), :), 'LineWidth', 2);
    %         hold on;
    %         no_color = no_color + 1;
    %     end
    %     no_color = 1;
    %     for ii = picked_snaps_contact
    %         plot(XY_tension_rotated(ii, :, 1), XY_tension_rotated(ii, :, 2), 'Color', ...
    %             color_group_contact(no_color+round(total_snapshots*unused_color_ratio/2), :), 'LineWidth', 2);
    %         hold on;
    %         no_color = no_color + 1;
    %     end
    %     XY_pillar_rotated(XY_pillar_rotated(:, 2) < 0.2e-4, :) = [];
    %     XY_pillar_rotated(XY_pillar_rotated(:, 2) > 1.3e-4, :) = [];
    %     XY_pillar_rotated(XY_pillar_rotated(:, 1) < 2.5e-4, :) = [];
    %     XY_pillar_rotated(XY_pillar_rotated(:, 1) > 5.5e-4, :) = [];
    %     viscircles(XY_pillar_rotated, R_pillar*ones(size(XY_pillar_rotated, 1), 1), 'Color', 'k', 'LineWidth', 0.2);
    %     axis equal; axis off
    %     xlim([3.3 4.7]*1e-4); ylim([5.5 9.5]*1e-5);
    %     set(gcf,'renderer','Painters');
    %     print('-depsc2','-tiff','-r100','-vector',[figure_mother_save_path, filesep, 'Long_fiber_chronophotography_new.eps']);
    % end
end





%% Load saved tension data, plot demo (SI Fig in the paper)
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

X_pillar = 0:C2C_pillar:1e-2; Y_pillar = 0:-C2C_pillar:-1e-2;
[X_pillar, Y_pillar] = meshgrid(X_pillar, Y_pillar);

theta = -35; % angle in degrees (clockwise)
Rotate_matrix = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; % rotation matrix

XY_pillar_rotated = [X_pillar(:), Y_pillar(:)]*Rotate_matrix;
% % figure('Position', [100, 100, 800, 800], 'Color', 'w');
% % hold on;
% % viscircles(XY_pillar_rotated, R_pillar*ones(size(XY_pillar_rotated, 1), 1), 'Color', 'r');
% % axis equal; axis off

figure_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Simulation\Tension analysis'];

parent_path = 'E:\Experimental Data (RAW)\Simulation\data';
sub1_path = dir(parent_path);

% Case: n90, no.1 (Long fiber)
sub1Path_i = 17; sub2Path_i = 6;

current_beadsNum = sub1_path(sub1Path_i).name;
beads_num = str2double(current_beadsNum(2:end));
disp(['Current beads number: ', num2str(beads_num)]);

L_0 = 0.6 * beads_num * 1e-6;
disp(['Current fiber length: ', num2str(L_0), ' um']);

sub2_path = dir(fullfile(parent_path, current_beadsNum));

current_sim = sub2_path(sub2Path_i).name;
sim_no = str2double(current_sim(11:end));
disp(['Current simulation number: ', num2str(sim_no)]);

fileinfo = dir(fullfile(parent_path, current_beadsNum, current_sim, 'output_data\*.vtk'));

picked_snap = 59;

snapshot = readVTK(fullfile(fileinfo(picked_snap).folder, fileinfo(picked_snap).name));
XY_current = snapshot.points(:, 1:2);
XY_current_rotated = XY_current*Rotate_matrix;

XY_current_in_lattice = mod(XY_current-C2C_pillar/2, C2C_pillar);
dist2origin = sqrt((XY_current_in_lattice(:,1)-C2C_pillar/2).^2 + (XY_current_in_lattice(:,2)-C2C_pillar/2).^2);
ContactBeads_index_current = dist2origin < R_pillar + 0.3 * 1e-6;

XY_tension_rotated = 0.5*(XY_current(1:end-1, :)+XY_current(2:end, :))*Rotate_matrix;

vec_im = diff(XY_current);
l_im = vecnorm(vec_im')';
vec_t = vec_im./l_im;
f_tension_magnitude_normalized = ks*(l_im - 2*R_fib)/f_bead; 


figure('Position', [100, 100, 800, 600], 'Color', 'w');

% plot pillar array
XY_pillar_rotated(XY_pillar_rotated(:, 2) < 0.2e-4, :) = [];
XY_pillar_rotated(XY_pillar_rotated(:, 2) > 1.3e-4, :) = [];
XY_pillar_rotated(XY_pillar_rotated(:, 1) < 2.5e-4, :) = [];
XY_pillar_rotated(XY_pillar_rotated(:, 1) > 5.5e-4, :) = [];
viscircles(XY_pillar_rotated, R_pillar*ones(size(XY_pillar_rotated, 1), 1), ...
    'Color', [0.3 0.3 0.3], 'LineWidth', 0.3, 'EnhanceVisibility', false);
axis equal; axis off
hold on;

% plot fiber
viscircles(XY_current_rotated(ContactBeads_index_current, :), R_fib*1.5, ...
    'Color', 'm', 'LineStyle', '-', 'LineWidth', 0.5); hold on;
hold on;

% plot tension
myColormap = cmocean('amp');
scatter(XY_tension_rotated(~ContactBeads_index_current(1:end-1), 1), ...
    XY_tension_rotated(~ContactBeads_index_current(1:end-1), 2), 50, ...
    f_tension_magnitude_normalized(~ContactBeads_index_current(1:end-1)), 'filled');
colormap(myColormap);
c = colorbar;
c.Label.String = '$\tilde{F}_{\rm{tension}}$';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 24;
c.TickLabelInterpreter = 'latex';
c.FontSize = 24;
clim([0 5]);
hold on;

% plot two fiber ends with different symbols
scatter(XY_current_rotated(1, 1), XY_current_rotated(1, 2), 100, 'k', 'filled', ...
    'd', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'k'); hold on;
scatter(XY_current_rotated(end, 1), XY_current_rotated(end, 2), 100, 'k', ...
    'filled', '>', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'k'); hold on;

% add a box around the figure
rectangle('Position', [3.4e-4, 5.6e-5, 0.7e-4, 3.8e-5], 'EdgeColor', 'k', 'LineWidth', 0.2); 
hold on;
xlim([3.4 4.1]*1e-4); ylim([5.6 9.4]*1e-5);

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',[figure_mother_save_path, filesep, 'theDEMO.eps']);