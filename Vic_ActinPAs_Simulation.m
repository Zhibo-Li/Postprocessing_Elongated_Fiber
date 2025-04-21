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
for sub1Path_i = 3:length(sub1_path)-1
    current_beadsNum = sub1_path(sub1Path_i).name;
    beads_num = str2double(current_beadsNum(2:end));

    L_0 = 0.6 * beads_num * 1e-6;

    sub2_path = dir(fullfile(parent_path, current_beadsNum));

    XY = cell(1500, length(sub2_path)-5); % 1500 is the maximum number of snapshots (roughly)
    CoM = cell(1500, length(sub2_path)-5); % 1500 is the maximum number of snapshots (roughly)
    Chi = cell(1500, length(sub2_path)-5); % 1500 is the maximum number of snapshots (roughly)
    aniso = cell(1500, length(sub2_path)-5); % 1500 is the maximum number of snapshots (roughly)
    L_ee_norm = cell(1500, length(sub2_path)-5); % 1500 is the maximum number of snapshots (roughly)
    Energy = cell(1500, length(sub2_path)-5); % 1500 is the maximum number of snapshots (roughly)
    f_tension = cell(1500, length(sub2_path)-5); % 1500 is the maximum number of snapshots (roughly)
    If_contact = cell(1500, length(sub2_path)-5); % 1500 is the maximum number of snapshots (roughly)


    for sub2Path_i = 6:length(sub2_path)
        current_sim = sub2_path(sub2Path_i).name;
        sim_no = str2double(current_sim(11:end));

        fileinfo = dir(fullfile(parent_path, current_beadsNum, current_sim, 'output_data\*.vtk'));

        for ii = 1:length(fileinfo)
            snapshot = readVTK(fullfile(fileinfo(ii).folder, fileinfo(ii).name));
            XY_current = snapshot.points(:, 1:2);

            XY_current_in_lattice = mod(XY_current-C2C_pillar/2, C2C_pillar);
            dist2origin = sqrt((XY_current_in_lattice(:,1)-C2C_pillar/2).^2 + (XY_current_in_lattice(:,2)-C2C_pillar/2).^2);
            if any(dist2origin < R_pillar + 0.6 * 1e-6)
                If_contact{ii, sub2Path_i-5} = 1;
            else
                If_contact{ii, sub2Path_i-5} = 0;
            end


            XY{ii, sub2Path_i-5} = snapshot.points(:, 1:2);

            vec_im = diff(XY_current);
            l_im = vecnorm(vec_im')';
            vec_t = vec_im./l_im;
            % f_tension_vector{ii, sub2Path_i-5} = -ks*(l_im - 2*R_fib).*vec_t;
            f_tension{ii, sub2Path_i-5} = ks*(l_im - 2*R_fib); % Notice the sign

            CoM_xy = mean(XY_current, 1);
            CoM{ii, sub2Path_i-5} = CoM_xy;
            Gyr = 1/size(XY_current,1) * [sum((XY_current(:, 1)-CoM_xy(1)).^2),  sum((XY_current(:, 1)-CoM_xy(1)) .* (XY_current(:, 2)-CoM_xy(2)));
                sum((XY_current(:, 2)-CoM_xy(2)) .* (XY_current(:, 1)-CoM_xy(1))), sum((XY_current(:, 2)-CoM_xy(2)).^2)];

            [eigenV,eigenD] = eig(Gyr);
            [d,ind] = sort(diag(eigenD));
            Ds = eigenD(ind,ind);
            Vs = eigenV(:,ind);
            Chi{ii, sub2Path_i-5} = atand(Vs(2,2)/Vs(1,2)); % orientation

            Lambda1 = eigenD(2,2); Lambda2 =  eigenD(1,1);
            aniso{ii, sub2Path_i-5} = 1 - 4*Lambda1*Lambda2/(Lambda1+Lambda2)^2; % sphericity

            L_ee_norm{ii, sub2Path_i-5} = sqrt((XY_current(1,1)-XY_current(end,1))^2 + (XY_current(1,2)-XY_current(end,2))^2) / L_0;

            [~,segl] = arclength(XY_current(:,1), XY_current(:,2), 'linear');
            [~,R,~] = curvature(XY_current); % calculate the curvature
            Energy{ii, sub2Path_i-5} = B / 2 * sum((1./R(2:end)).^2 .* segl, 'omitnan');

            % for foo = 1:50
            %     plot(XY_current(foo,1), XY_current(foo,2), 'm', 'Marker','o', 'MarkerSize', 10); hold on;
            %     axis equal; axis off;
            % end

            % if mod(ii, 10) == 0
            %     plot(XY_current(:,1), XY_current(:,2), '-m'); hold on;
            %     axis equal; axis off;
            % end

        end
    end

    % remove the empty cell
    XY = XY(~cellfun('isempty', XY));
    f_tension = f_tension(~cellfun('isempty', f_tension));
    CoM = CoM(~cellfun('isempty', CoM));
    Chi = Chi(~cellfun('isempty', Chi));
    aniso = aniso(~cellfun('isempty', aniso));
    L_ee_norm = L_ee_norm(~cellfun('isempty', L_ee_norm));
    Energy = Energy(~cellfun('isempty', Energy));
    If_contact = logical(cell2mat(If_contact(~cellfun('isempty', If_contact))));

    if ~exist("data_mother_save_path", 'dir')
        mkdir(data_mother_save_path);
    end

    save([data_mother_save_path, filesep, 'FiberLength=', num2str(L_0*1e6), 'um.mat'], ...
        'XY', 'f_tension', 'CoM', 'Chi', 'aniso', 'L_ee_norm', 'Energy', 'If_contact', 'L_0');

end





%% Load saved data and plot the results
clear; close all; clc;

% Geometrical properties of the array
R_pillar = 8.3 * 1e-6; % fiber radius, unit: m
C2C_pillar = 30 * 1e-6; % m (delta_x & delta_y)
delta_t = 0.1; % s

Xedges = 0:0.02:1; Yedges = 0:0.02:1; % for the histogram


data_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Simulation\Simulations_LubricationModel'];
saved_data = dir(data_mother_save_path);

figure_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Simulation\Velocity analysis'];

L_0_and_num_switch_head_percentage = zeros(length(saved_data)-2, 2);
L_0_and_num_contact_percentage = zeros(length(saved_data)-2, 2);
L_0_and_num_only_ends_contact_percentage = zeros(length(saved_data)-2, 2);
for ii = 3:length(saved_data)

    load(fullfile(saved_data(ii).folder, saved_data(ii).name));

    L_0 = extractBetween(saved_data(ii).name, 'FiberLength=', 'um.mat'); % unit: um
    bead_num = size(XY{1, 1}, 1);

    % Plot the velocity of the fiber/fiber ends vs. time/position
    XY_all = cell2mat(XY);
    last_bead_pos = XY_all(1:bead_num:end, :);
    first_bead_pos = XY_all(bead_num:bead_num:end, :);
    CoM_all = cell2mat(CoM);
    trajectory_cut_index = [0; find(abs(diff(CoM_all(:,1))) > 2e-3); length(CoM_all)];
    for jj = 1:10:length(trajectory_cut_index)-1
        XY_current = XY_all(trajectory_cut_index(jj)+1:trajectory_cut_index(jj+1), :);

        CoM_current = CoM_all(trajectory_cut_index(jj)+1:trajectory_cut_index(jj+1), :);
        CoM_velocity = diff(CoM_current) / delta_t;
        % convert the velocity to the local coordinate system (x is the main flow direction)
        CoM_velocity_rotated = CoM_velocity * [cosd(35), -sind(35); sind(35), cosd(35)];
        CoM_velocity_rotated_norm = sqrt(sum(CoM_velocity_rotated.^2, 2));
        figure('Position', [100, 100, 2400, 200], 'Color', 'w');
        plot(CoM_current(1:end-1, 1), CoM_velocity_rotated_norm); hold on;
        % plot the position of the pillars
        pillar_positions = 0:C2C_pillar:max(CoM_current(:, 1));
        arrayfun(@(x) line([x x], [0 2e-4], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 0.2), pillar_positions);
        xlim([0 3e-3]);ylim([0 2e-4]);
        ylabel('Velocity of the CoM (um/s)');
        xlabel('Position (um)');
        % saveas(gcf, [figure_mother_save_path, filesep, 'Velocity of the CoM, L = ', L_0{1}, ' um - ', num2str(jj), '.png']);
        
        % Fourier analysis of CoM velocity
        L = length(CoM_velocity_rotated_norm);
        Y = fft(CoM_velocity_rotated_norm);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = (0:(L/2))/L/delta_t;
        
        figure('Position', [100, 100, 1200, 400], 'Color', 'w');
        plot(f, P1);
        title('Single-Sided Amplitude Spectrum of CoM Velocity');
        xlabel('Frequency (Hz)');
        ylabel('|P1(f)|');
        xlim([0 1/delta_t/2]);
        grid on;
        % saveas(gcf, [figure_mother_save_path, filesep, 'Fourier Analysis of CoM Velocity, L = ', L_0{1}, ' um - ', num2str(jj), '.png']);

        first_bead_pos_current = first_bead_pos(trajectory_cut_index(jj)+1:trajectory_cut_index(jj+1), :);
        last_bead_pos_current = last_bead_pos(trajectory_cut_index(jj)+1:trajectory_cut_index(jj+1), :);
        first_bead_pos_velocity = diff(first_bead_pos_current) / delta_t;
        last_bead_pos_velocity = diff(last_bead_pos_current) / delta_t;
        % convert the velocity to the local coordinate system (x is the main flow direction)
        first_bead_pos_velocity_rotated = first_bead_pos_velocity * [cosd(35), -sind(35); sind(35), cosd(35)];
        last_bead_pos_velocity_rotated = last_bead_pos_velocity * [cosd(35), -sind(35); sind(35), cosd(35)];
        first_bead_pos_velocity_rotated_norm = sqrt(sum(first_bead_pos_velocity_rotated.^2, 2));
        last_bead_pos_velocity_rotated_norm = sqrt(sum(last_bead_pos_velocity_rotated.^2, 2));
        figure('Position', [100, 100, 2400, 200], 'Color', 'w');
        % plot the position of the pillars
        arrayfun(@(x) line([x x], [0 2e-4], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 0.2), pillar_positions); hold on;
        yyaxis left
        plot(first_bead_pos_current(1:end-1, 1), first_bead_pos_velocity_rotated_norm, 'b'); hold on;
        ylabel('Velocity of the first bead (um/s)'); ylim([0 2e-4]);
        yyaxis right
        plot(last_bead_pos_current(1:end-1, 1), last_bead_pos_velocity_rotated_norm, 'g'); hold on;
        ylabel('Velocity of the last bead (um/s)'); ylim([0 2e-4]);
        xlim([0 3e-3]); 
        xlabel('Position (um)');
        % saveas(gcf, [figure_mother_save_path, filesep, 'Velocity of the fiber ends, L = ', L_0{1}, ' um - ', num2str(jj), '.png']);

    end



    % % Plot the normalized fiber tension in the unit cell
    % f_tension_XY = cellfun(@(x) 0.5*(x(1:end-1, :)+x(2:end, :)), XY, 'UniformOutput', false);
    % f_tension_all = cell2mat(f_tension);
    % f_tension_all_XY = cell2mat(f_tension_XY);
    % f_tension_XY_in_lattice = mod(f_tension_all_XY-C2C_pillar/2, C2C_pillar);
    % f_tension_XY_in_lattice_normalized = f_tension_XY_in_lattice ./ [C2C_pillar, C2C_pillar];
    % f_tension_in_grid_norm_mean = zeros(length(Xedges)-1, length(Yedges)-1);
    % for grid_x = 1:length(Xedges)-1
    %     for grid_y = 1:length(Yedges)-1
    %         If_fiber_in_grid = f_tension_XY_in_lattice_normalized(:,1) >= Xedges(grid_x) & ...
    %             f_tension_XY_in_lattice_normalized(:,1) < Xedges(grid_x+1) & ...
    %             f_tension_XY_in_lattice_normalized(:,2) >= Yedges(grid_y) & ...
    %             f_tension_XY_in_lattice_normalized(:,2) < Yedges(grid_y+1);
    %         f_tension_in_grid = f_tension_all(If_fiber_in_grid);
    %         f_tension_in_grid_norm = sqrt(sum(f_tension_in_grid.^2, 2));
    %         f_tension_in_grid_norm_mean(grid_x, grid_y) = mean(f_tension_in_grid_norm);
    %     end
    % end
    % figure('Position', [100, 100, 800, 800], 'Color', 'w');
    % [XX, YY] = meshgrid(movmean(Xedges, 2, "Endpoints","discard"), movmean(Yedges, 2, "Endpoints","discard"));
    % % Notice the transpose of the matrix!!!
    % pcolor(XX, YY, f_tension_in_grid_norm_mean', "EdgeColor", "none");
    % axis equal; axis off; grid off
    % viscircles([0.5 0.5], R_pillar/C2C_pillar, 'Color', 'r', 'LineStyle', '--');
    % xlim([0 1]); ylim([0 1]);
    % cmocean('thermal'); colorbar; clim([0 4e-9]);
    % title('Mean of the fiber tension norm, L = ' + string(L_0) + ' um');
    % saveas(gcf, [figure_mother_save_path, filesep, 'Mean of the fiber tension norm, L = ', L_0{1}, ' um - All.png']);



    % % Plot the fiber and CoM in the unit cell (when one ends of the fibers contact pillars)
    % contact_threshold = 0.3 * 1e-6;
    % If_contact_re_calcu = zeros(length(XY), 1);
    % If_rest_beads_contact = zeros(length(XY), 1);
    % contact_bead_pos = cell(length(XY), 1);
    % for fiber_i = 1:length(XY)
    %     current_fiber_XY = XY{fiber_i};
    %     current_fiber_XY_in_lattice = mod(current_fiber_XY-C2C_pillar/2, C2C_pillar);
    %     if any(sqrt(sum((current_fiber_XY_in_lattice - C2C_pillar/2).^2, 2)) < R_pillar + contact_threshold)
    %         If_contact_re_calcu(fiber_i) = 1;
    %         current_fiber_NO_ends_XY = current_fiber_XY(2:end-1, :);
    %         current_fiber_NO_ends_XY_in_lattice = mod(current_fiber_NO_ends_XY-C2C_pillar/2, C2C_pillar);
    %         if any(sqrt(sum((current_fiber_NO_ends_XY_in_lattice - C2C_pillar/2).^2, 2)) < R_pillar + contact_threshold)
    %             If_rest_beads_contact(fiber_i) = 1;
    %         else
    %             If_rest_beads_contact(fiber_i) = 0;
    %             last_bead_pos = current_fiber_XY(1, :); last_bead_pos_in_lattice = mod(last_bead_pos-C2C_pillar/2, C2C_pillar);
    %             first_bead_pos = current_fiber_XY(end, :); first_bead_pos_in_lattice = mod(first_bead_pos-C2C_pillar/2, C2C_pillar);
    %             if sum((last_bead_pos_in_lattice - C2C_pillar/2).^2) < sum((first_bead_pos_in_lattice - C2C_pillar/2).^2)
    %                 contact_bead_pos{fiber_i} = last_bead_pos;
    %             else
    %                 contact_bead_pos{fiber_i} = first_bead_pos;
    %             end
    %         end
    %     else
    %         If_contact_re_calcu(fiber_i) = 0;
    %     end
    % end
    % XY_only_ends_contact = XY(If_contact_re_calcu & ~If_rest_beads_contact); % fibers that only one end contacts the pillar
    % XY_only_ends_contact_in_lattice = mod(cell2mat(XY_only_ends_contact)-C2C_pillar/2, C2C_pillar);
    % XY_only_ends_contact_in_lattice_normalized = XY_only_ends_contact_in_lattice ./ [C2C_pillar, C2C_pillar]; % normalized to [0, 1]
    % CoM_only_ends_contact = CoM(If_contact_re_calcu & ~If_rest_beads_contact); % CoM of fibers that only one end contacts the pillar
    % % CoM_only_ends_contact_in_lattice = mod(cell2mat(CoM_only_ends_contact)-C2C_pillar/2, C2C_pillar);
    % % CoM_only_ends_contact_in_lattice_normalized = CoM_only_ends_contact_in_lattice ./ [C2C_pillar, C2C_pillar]; % normalized to [0, 1]
    % contact_bead_pos_all = contact_bead_pos(If_contact_re_calcu & ~If_rest_beads_contact); % contact bead position of fibers that only one end contacts the pillar
    % contact_bead_pos_all_in_lattice = mod(cell2mat(contact_bead_pos_all)-C2C_pillar/2, C2C_pillar);
    % contact_bead_pos_all_in_lattice_normalized = contact_bead_pos_all_in_lattice ./ [C2C_pillar, C2C_pillar]; % normalized to [0, 1]
    % CoM_only_ends_contact_rotated = cell2mat(CoM_only_ends_contact) * [cosd(-35), -sind(-35); sind(-35), cosd(-35)];
    % contact_bead_pos_all_in_lattice_normalized_rotated = (contact_bead_pos_all_in_lattice_normalized - [0, 1]) * [cosd(-35), -sind(-35); sind(-35), cosd(-35)]; % Rotate based on the upper left corner
    % contact_bead_pos_all_rotated = cell2mat(contact_bead_pos_all) * [cosd(-35), -sind(-35); sind(-35), cosd(-35)];
    % fiber_leave_pillar_index = CoM_only_ends_contact_rotated(:,1) > contact_bead_pos_all_rotated(:,1);
    % fiber_above_leave_pillar_index = contact_bead_pos_all_in_lattice_normalized_rotated(:, 2) > 0;
    % % figure('Position', [100, 100, 800, 800], 'Color', 'w');
    % % for fiber_i = 1:length(XY_only_ends_contact)
    % %     plot(mod(XY_only_ends_contact{fiber_i}(:,1)-C2C_pillar/2, C2C_pillar) / ...
    % %         C2C_pillar, mod(XY_only_ends_contact{fiber_i}(:,2)-C2C_pillar/2, C2C_pillar) / ...
    % %         C2C_pillar, 'm', 'LineStyle', 'none', 'Marker', '*', 'MarkerSize', 3); hold on;
    % %     plot(mod(CoM_only_ends_contact{fiber_i}(1)-C2C_pillar/2, C2C_pillar) / ...
    % %         C2C_pillar, mod(CoM_only_ends_contact{fiber_i}(2)-C2C_pillar/2, C2C_pillar) / ...
    % %         C2C_pillar, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
    % % end
    % % axis equal; axis off;
    % % viscircles([0.5 0.5], R_pillar/C2C_pillar, 'Color', 'r', 'LineStyle', '--');
    % % xlim([0 1]); ylim([0 1]);
    % % title('Fibers (only one end contact) and their CoM, L = ' + string(L_0) + ' um, All');
    % % saveas(gcf, [figure_mother_save_path, filesep, 'Fibers (only one end contact) and their CoM, L = ', L_0{1}, ' um - All.png']);

    % XY_only_ends_contact_and_fiber_above_leave_pillar = XY_only_ends_contact(fiber_leave_pillar_index & fiber_above_leave_pillar_index);
    % CoM_only_ends_contact_and_fiber_above_leave_pillar = CoM_only_ends_contact(fiber_leave_pillar_index & fiber_above_leave_pillar_index);
    % XY_only_ends_contact_and_fiber_below_leave_pillar = XY_only_ends_contact(fiber_leave_pillar_index & ~fiber_above_leave_pillar_index);
    % CoM_only_ends_contact_and_fiber_below_leave_pillar = CoM_only_ends_contact(fiber_leave_pillar_index & ~fiber_above_leave_pillar_index);
    % figure('Position', [100, 100, 800, 800], 'Color', 'w');
    % for fiber_i = 1:length(XY_only_ends_contact_and_fiber_above_leave_pillar)
    %     plot(mod(XY_only_ends_contact_and_fiber_above_leave_pillar{fiber_i}(:,1)-C2C_pillar/2, C2C_pillar) / ...
    %         C2C_pillar, mod(XY_only_ends_contact_and_fiber_above_leave_pillar{fiber_i}(:,2)-C2C_pillar/2, C2C_pillar) / ...
    %         C2C_pillar, 'm', 'LineStyle', 'none', 'Marker', '*', 'MarkerSize', 3); hold on;
    %     plot(mod(CoM_only_ends_contact_and_fiber_above_leave_pillar{fiber_i}(1)-C2C_pillar/2, C2C_pillar) / ...
    %         C2C_pillar, mod(CoM_only_ends_contact_and_fiber_above_leave_pillar{fiber_i}(2)-C2C_pillar/2, C2C_pillar) / ...
    %         C2C_pillar, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); hold on;
    % end
    % for fiber_i = 1:length(XY_only_ends_contact_and_fiber_below_leave_pillar)
    %     plot(mod(XY_only_ends_contact_and_fiber_below_leave_pillar{fiber_i}(:,1)-C2C_pillar/2, C2C_pillar) / ...
    %         C2C_pillar, mod(XY_only_ends_contact_and_fiber_below_leave_pillar{fiber_i}(:,2)-C2C_pillar/2, C2C_pillar) / ...
    %         C2C_pillar, 'c', 'LineStyle', 'none', 'Marker', '*', 'MarkerSize', 3); hold on;
    %     plot(mod(CoM_only_ends_contact_and_fiber_below_leave_pillar{fiber_i}(1)-C2C_pillar/2, C2C_pillar) / ...
    %         C2C_pillar, mod(CoM_only_ends_contact_and_fiber_below_leave_pillar{fiber_i}(2)-C2C_pillar/2, C2C_pillar) / ...
    %         C2C_pillar, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); hold on;
    % end
    % axis equal; axis off;
    % viscircles([0.5 0.5], R_pillar/C2C_pillar, 'Color', 'k', 'LineStyle', '--');
    % xlim([0 1]); ylim([0 1]);
    % rectangle('Position', [0, 0, 1, 1], 'EdgeColor', 'k', 'LineWidth', 0.5);
    % title('Fibers (only one end contact and fiber leaves pillar) and their CoM, L = ' + string(L_0) + ' um, All');
    % saveas(gcf, [figure_mother_save_path, filesep, 'Fibers (only one end contact and fiber leaves pillar) and their CoM, L = ', L_0{1}, ' um - All.png']);



    % % Plot the fiber tension along its centerline in the unit cell
    % contact_threshold = 0.3 * 1e-6;
    % If_contact_re_calcu = zeros(length(XY), 1);
    % If_rest_beads_contact = zeros(length(XY), 1);
    % contact_bead_pos = cell(length(XY), 1);
    % for fiber_i = 1:length(XY)
    %     current_fiber_XY = XY{fiber_i};
    %     current_fiber_XY_in_lattice = mod(current_fiber_XY-C2C_pillar/2, C2C_pillar);
    %     if any(sqrt(sum((current_fiber_XY_in_lattice - C2C_pillar/2).^2, 2)) < R_pillar + contact_threshold)
    %         If_contact_re_calcu(fiber_i) = 1;
    %         current_fiber_NO_ends_XY = current_fiber_XY(2:end-1, :);
    %         current_fiber_NO_ends_XY_in_lattice = mod(current_fiber_NO_ends_XY-C2C_pillar/2, C2C_pillar);
    %         if any(sqrt(sum((current_fiber_NO_ends_XY_in_lattice - C2C_pillar/2).^2, 2)) < R_pillar + contact_threshold)
    %             If_rest_beads_contact(fiber_i) = 1;
    %         else
    %             If_rest_beads_contact(fiber_i) = 0;
    %             last_bead_pos = current_fiber_XY(1, :); last_bead_pos_in_lattice = mod(last_bead_pos-C2C_pillar/2, C2C_pillar);
    %             first_bead_pos = current_fiber_XY(end, :); first_bead_pos_in_lattice = mod(first_bead_pos-C2C_pillar/2, C2C_pillar);
    %             if sum((last_bead_pos_in_lattice - C2C_pillar/2).^2) < sum((first_bead_pos_in_lattice - C2C_pillar/2).^2)
    %                 contact_bead_pos{fiber_i} = last_bead_pos;
    %             else
    %                 contact_bead_pos{fiber_i} = first_bead_pos;
    %             end
    %         end
    %     else
    %         If_contact_re_calcu(fiber_i) = 0;
    %     end
    % end
    % XY_only_ends_contact = XY(If_contact_re_calcu & ~If_rest_beads_contact); % fibers that only one end contacts the pillar
    % % XY_only_ends_contact_in_lattice = mod(cell2mat(XY_only_ends_contact)-C2C_pillar/2, C2C_pillar);
    % % XY_only_ends_contact_in_lattice_normalized = XY_only_ends_contact_in_lattice ./ [C2C_pillar, C2C_pillar]; % normalized to [0, 1]
    % CoM_only_ends_contact = CoM(If_contact_re_calcu & ~If_rest_beads_contact); % CoM of fibers that only one end contacts the pillar
    % % CoM_only_ends_contact_in_lattice = mod(cell2mat(CoM_only_ends_contact)-C2C_pillar/2, C2C_pillar);
    % % CoM_only_ends_contact_in_lattice_normalized = CoM_only_ends_contact_in_lattice ./ [C2C_pillar, C2C_pillar]; % normalized to [0, 1]
    % contact_bead_pos_all = contact_bead_pos(If_contact_re_calcu & ~If_rest_beads_contact); % contact bead position of fibers that only one end contacts the pillar
    % % contact_bead_pos_all_in_lattice = mod(cell2mat(contact_bead_pos_all)-C2C_pillar/2, C2C_pillar);
    % % contact_bead_pos_all_in_lattice_normalized = contact_bead_pos_all_in_lattice ./ [C2C_pillar, C2C_pillar]; % normalized to [0, 1]
    % f_tension_only_ends_contact = f_tension(If_contact_re_calcu & ~If_rest_beads_contact);

    % CoM_only_ends_contact_rotated = cell2mat(CoM_only_ends_contact) * [cosd(-35), -sind(-35); sind(-35), cosd(-35)];
    % contact_bead_pos_all_rotated = cell2mat(contact_bead_pos_all) * [cosd(-35), -sind(-35); sind(-35), cosd(-35)];
    % fiber_leave_pillar_index = CoM_only_ends_contact_rotated(:,1) > contact_bead_pos_all_rotated(:,1);
    % XY_only_ends_contact_and_fiber_leave_pillar = XY_only_ends_contact(fiber_leave_pillar_index);
    % CoM_only_ends_contact_and_fiber_leave_pillar = CoM_only_ends_contact(fiber_leave_pillar_index);
    % f_tension_only_ends_contact_and_fiber_leave_pillar = f_tension_only_ends_contact(fiber_leave_pillar_index);
    % figure('Position', [100, 100, 800, 800], 'Color', 'w');
    % for fiber_i = 1:length(XY_only_ends_contact_and_fiber_leave_pillar)
    %     scatter(mod(movmean(XY_only_ends_contact_and_fiber_leave_pillar{fiber_i}(:,1), 2, "Endpoints", "discard")-C2C_pillar/2, C2C_pillar) / ...
    %         C2C_pillar, mod(movmean(XY_only_ends_contact_and_fiber_leave_pillar{fiber_i}(:,2), 2, "Endpoints", "discard")-C2C_pillar/2, C2C_pillar) / ...
    %         C2C_pillar, 24, sqrt(sum(f_tension_only_ends_contact_and_fiber_leave_pillar{fiber_i}.^2, 2))); hold on;
    % end
    % axis equal; axis off;
    % viscircles([0.5 0.5], R_pillar/C2C_pillar, 'Color', 'r', 'LineStyle', '--');
    % xlim([0 1]); ylim([0 1]);
    % title('Fiber tension along its centerline, L = ' + string(L_0) + ' um, All');

    % figure('Position', [100, 100, 800, 800], 'Color', 'w');
    % for fiber_i = 1:length(XY_only_ends_contact_and_fiber_leave_pillar)
    %     plot(sqrt(sum(f_tension_only_ends_contact_and_fiber_leave_pillar{fiber_i}.^2, 2))); hold on;
    % end
    % title('Fiber tension along its centerline, L = ' + string(L_0) + ' um, All');



    % % Plot probability distribution function of the polymer position
    % figure('Position', [100, 100, 800, 800], 'Color', 'w');
    % XY_all = cell2mat(XY);
    % fiber_in_lattice = mod(XY_all-C2C_pillar/2, C2C_pillar);
    % fiber_in_lattice_normalized = fiber_in_lattice ./ [C2C_pillar, C2C_pillar];
    % h = histogram2(fiber_in_lattice_normalized(:,1), fiber_in_lattice_normalized(:,2), Xedges, ...
    %     Yedges, 'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
    % axis equal; axis off
    % cmocean('ice'); % clim([0 10]); % colorbar;
    % viscircles([0.5 0.5], R_pillar/C2C_pillar, 'Color', 'r', 'LineStyle', '--');
    % xlim([0 1]); ylim([0 1]);
    % title('PDF of the polymer position, L = ' + string(L_0) + ' um, All');
    % saveas(gcf, [figure_mother_save_path, filesep, 'PDF of the polymer position, L = ', L_0{1}, ' um - All.png']);



    % % Plot probability distribution function of the polymer position (when fibers contact pillars)
    % figure('Position', [100, 100, 800, 800], 'Color', 'w');
    % XY_all = cell2mat(XY(If_contact));
    % fiber_in_lattice = mod(XY_all-C2C_pillar/2, C2C_pillar);
    % fiber_in_lattice_normalized = fiber_in_lattice ./ [C2C_pillar, C2C_pillar];
    % h = histogram2(fiber_in_lattice_normalized(:,1), fiber_in_lattice_normalized(:,2), Xedges, ...
    %     Yedges, 'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
    % axis equal; axis off
    % cmocean('ice'); % clim([0 10]); % colorbar;
    % viscircles([0.5 0.5], R_pillar/C2C_pillar, 'Color', 'r', 'LineStyle', '--');
    % xlim([0 1]); ylim([0 1]);
    % title('PDF of the polymer position, L = ' + string(L_0) + ' um, Contact');
    % saveas(gcf, [figure_mother_save_path, filesep, 'PDF of the polymer position, L = ', L_0{1}, ' um - Contact.png']);



    % % Plot probability distribution function of the polymer position (when fibers don't contact pillars)
    % figure('Position', [100, 100, 800, 800], 'Color', 'w');
    % XY_all = cell2mat(XY(~If_contact));
    % fiber_in_lattice = mod(XY_all-C2C_pillar/2, C2C_pillar);
    % fiber_in_lattice_normalized = fiber_in_lattice ./ [C2C_pillar, C2C_pillar];
    % h = histogram2(fiber_in_lattice_normalized(:,1), fiber_in_lattice_normalized(:,2), Xedges, ...
    %     Yedges, 'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
    % axis equal; axis off
    % cmocean('ice'); % clim([0 10]); % colorbar;
    % viscircles([0.5 0.5], R_pillar/C2C_pillar, 'Color', 'r', 'LineStyle', '--');
    % xlim([0 1]); ylim([0 1]);
    % title('PDF of the polymer position, L = ' + string(L_0) + ' um, No contact');
    % saveas(gcf, [figure_mother_save_path, filesep, 'PDF of the polymer position, L = ', L_0{1}, ' um - No contact.png']);



    % % plot probability distribution function of the polymer COM
    % figure('Position', [100, 100, 800, 800], 'Color', 'w');
    % CoM_in_lattice = mod(cell2mat(CoM)-C2C_pillar/2, C2C_pillar);
    % CoM_in_lattice_normalized = CoM_in_lattice ./ [C2C_pillar, C2C_pillar];
    % h = histogram2(CoM_in_lattice_normalized(:,1), CoM_in_lattice_normalized(:,2), Xedges, ...
    %     Yedges, 'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
    % axis equal; axis off
    % cmocean('thermal'); % clim([0 10]); % colorbar;
    % viscircles([0.5 0.5], R_pillar/C2C_pillar, 'Color', 'r', 'LineStyle', '--');
    % xlim([0 1]); ylim([0 1]);
    % title('PDF of the polymer CoM, L = ' + string(L_0) + ' um, All');
    % saveas(gcf, [figure_mother_save_path, filesep, 'PDF of the polymer CoM, L = ', L_0{1}, ' um - All.png']);



    % % plot probability distribution function of the polymer COM (when fibers contact pillars)
    % figure('Position', [100, 100, 800, 800], 'Color', 'w');
    % CoM_Contact = CoM(If_contact);
    % CoM_in_lattice = mod(cell2mat(CoM_Contact)-C2C_pillar/2, C2C_pillar);
    % CoM_in_lattice_normalized = CoM_in_lattice ./ [C2C_pillar, C2C_pillar];
    % h = histogram2(CoM_in_lattice_normalized(:,1), CoM_in_lattice_normalized(:,2), Xedges, ...
    %     Yedges, 'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
    % axis equal; axis off
    % cmocean('thermal'); % clim([0 10]); % colorbar;
    % viscircles([0.5 0.5], R_pillar/C2C_pillar, 'Color', 'r', 'LineStyle', '--');
    % xlim([0 1]); ylim([0 1]);
    % title('PDF of the polymer CoM, L = ' + string(L_0) + ' um, Contact');
    % saveas(gcf, [figure_mother_save_path, filesep, 'PDF of the polymer CoM, L = ', L_0{1}, ' um - Contact.png']);



    % % plot probability distribution function of the polymer COM (when fibers don't contact pillars)
    % figure('Position', [100, 100, 800, 800], 'Color', 'w');
    % CoM_NoContact = CoM(~If_contact);
    % CoM_in_lattice = mod(cell2mat(CoM_NoContact)-C2C_pillar/2, C2C_pillar);
    % CoM_in_lattice_normalized = CoM_in_lattice ./ [C2C_pillar, C2C_pillar];
    % h = histogram2(CoM_in_lattice_normalized(:,1), CoM_in_lattice_normalized(:,2), Xedges, ...
    %     Yedges, 'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
    % axis equal; axis off
    % cmocean('thermal'); % clim([0 10]); % colorbar;
    % viscircles([0.5 0.5], R_pillar/C2C_pillar, 'Color', 'r', 'LineStyle', '--');
    % xlim([0 1]); ylim([0 1]);
    % title('PDF of the polymer CoM, L = ' + string(L_0) + ' um, No contact');
    % saveas(gcf, [figure_mother_save_path, filesep, 'PDF of the polymer CoM, L = ', L_0{1}, ' um - No contact.png']);



    % % L_0 and the ratio of the number of fibers that switch the head direction to the total number of fibers
    % XY_all = cell2mat(XY);
    % % rotate XY 35 drgrees counterclockwise
    % XY_all_rotated = XY_all * [cosd(-35), -sind(-35); sind(-35), cosd(-35)];
    % last_bead_pos = XY_all_rotated(1:bead_num:end, :);
    % first_bead_pos = XY_all_rotated(bead_num:bead_num:end, :);
    % % if fisrt_bead_pos is on the right side of the last_bead_pos, if_fiber_head_foward = 1; otherwise, if_fiber_head_foward = 0
    % if_fiber_head_foward = first_bead_pos(:, 1) > last_bead_pos(:, 1);
    % % figure('Position', [100, 100, 800, 800], 'Color', 'w');
    % % plot(if_fiber_head_foward, 'o', 'MarkerSize', 1, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'm');
    % % title('Fiber head direction, L = ' + string(L_0) + ' um');
    % % ylim([-0.5 1.5]);
    % % the number of fibers that switch the head direction
    % num_switch_head = sum(abs(diff(if_fiber_head_foward)) == 1);
    % % the ratio of the number of fibers that switch the head direction to the total number of fibers
    % num_switch_head_percentage = num_switch_head / length(if_fiber_head_foward);
    % % L_0 and num_switch_head_percentage into a matrix
    % L_0_and_num_switch_head_percentage(ii-2, :) = [str2double(L_0{1}), num_switch_head_percentage];



    % % % L_0 and the percentage of fibers that contact pillars
    % contact_threshold = 0.3 * 1e-6;
    % If_contact_re_calcu = zeros(length(XY), 1);
    % If_rest_beads_contact = zeros(length(XY), 1);
    % contact_bead_pos = cell(length(XY), 1);
    % for fiber_i = 1:length(XY)
    %     current_fiber_XY = XY{fiber_i};
    %     current_fiber_XY_in_lattice = mod(current_fiber_XY-C2C_pillar/2, C2C_pillar);
    %     if any(sqrt(sum((current_fiber_XY_in_lattice - C2C_pillar/2).^2, 2)) < R_pillar + contact_threshold)
    %         If_contact_re_calcu(fiber_i) = 1;
    %         current_fiber_NO_ends_XY = current_fiber_XY(2:end-1, :);
    %         current_fiber_NO_ends_XY_in_lattice = mod(current_fiber_NO_ends_XY-C2C_pillar/2, C2C_pillar);
    %         if any(sqrt(sum((current_fiber_NO_ends_XY_in_lattice - C2C_pillar/2).^2, 2)) < R_pillar + contact_threshold)
    %             If_rest_beads_contact(fiber_i) = 1;
    %         else
    %             If_rest_beads_contact(fiber_i) = 0;
    %             last_bead_pos = current_fiber_XY(1, :); last_bead_pos_in_lattice = mod(last_bead_pos-C2C_pillar/2, C2C_pillar);
    %             first_bead_pos = current_fiber_XY(end, :); first_bead_pos_in_lattice = mod(first_bead_pos-C2C_pillar/2, C2C_pillar);
    %             if sum((last_bead_pos_in_lattice - C2C_pillar/2).^2) < sum((first_bead_pos_in_lattice - C2C_pillar/2).^2)
    %                 contact_bead_pos{fiber_i} = last_bead_pos;
    %             else
    %                 contact_bead_pos{fiber_i} = first_bead_pos;
    %             end
    %         end
    %     else
    %         If_contact_re_calcu(fiber_i) = 0;
    %     end
    % end
    % % L_0 and percentage of fibers that contact pillars
    % L_0_and_num_contact_percentage(ii-2, :) = [str2double(L_0{1}), sum(If_contact_re_calcu) / length(If_contact_re_calcu)];
    % % L_0 and percentage of fibers that only one end contacts the pillar
    % L_0_and_num_only_ends_contact_percentage(ii-2, :) = [str2double(L_0{1}), sum(If_contact_re_calcu & ~If_rest_beads_contact) / length(If_contact_re_calcu)];


    clearvars XY CoM Chi aniso L_ee_norm Energy If_contact
    close all

end

% figure('Position', [100, 100, 800, 600], 'Color', 'w');
% plot(L_0_and_num_switch_head_percentage(:, 1), L_0_and_num_switch_head_percentage(:, 2)*100, ...
%     'o', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'm');
% ytickformat('percentage');
% xlabel('$L_0$ ($\mu$m)', 'Interpreter', 'latex', 'FontSize', 20);
% ylabel('Fiber switch-head percentage', 'Interpreter', 'latex', 'FontSize', 20);
% title('Fiber switch-head percentage vs. $L_0$', 'Interpreter', 'latex', 'FontSize', 20);
% set(gca, 'FontSize', 20);
% grid on;
% saveas(gcf, [figure_mother_save_path, filesep, 'Fiber switch-head percentage vs L_0.png']);

% figure('Position', [100, 100, 800, 600], 'Color', 'w');
% yyaxis left;
% plot(L_0_and_num_contact_percentage(:, 1), L_0_and_num_contact_percentage(:, 2)*100, ...
%     '^', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); hold on;
% ylabel('Fiber contact', 'Interpreter', 'latex', 'FontSize', 20);
% ytickformat('percentage');
% yyaxis right;
% plot(L_0_and_num_only_ends_contact_percentage(:, 1), L_0_and_num_only_ends_contact_percentage(:, 2)*100, ...
%     's', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
% ylabel('Fiber only the ends contact', 'Interpreter', 'latex', 'FontSize', 20);
% ytickformat('percentage');
% xlabel('$L_0$ ($\mu$m)', 'Interpreter', 'latex', 'FontSize', 20);
% title('Fiber contact percentage vs. $L_0$', 'Interpreter', 'latex', 'FontSize', 20);
% set(gca, 'FontSize', 20);
% grid on;
% saveas(gcf, [figure_mother_save_path, filesep, 'Fiber contact percentage vs L_0.png']);
