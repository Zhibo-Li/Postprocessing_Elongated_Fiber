%%% actin in PAs: DLD check
%
% lambda = 30 um, exp date: 20231030, 20240415, 20241031
% lambda = 45 um, exp date: 20240522, 20241111
%
% X_base = 6000 um: exp date: 20241031, 20241111
% X_base = 9000 um: exp date: 20231030, 20240415, 20240522
%
% flowfocusing_bandwidth = 100 um: exp date: 20241031, 20241111
% flowfocusing_bandwidth = 150 um: exp date: 20231030, 20240415, 20240522
% 
% NOTE: line59 uses exp_date_idx = [2,4,6,7] for X_base=6000um + Long fiber
%       experiments.

clear; close all; clc;

image_height = 204.8; % the height of the image, unit: um
pixel_size = 0.1; % the magnification of the objective (um/pixel)

set(0, 'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

data_path = 'Z:\Processing & Results\Actin Filaments in Porous Media\Actin-DLD';

subfolders = dir(data_path);
subfolders = subfolders([subfolders.isdir] & ~ismember({subfolders.name}, {'.', '..'}));
subfolders = subfolders(~cellfun(@isempty, regexp({subfolders.name}, '\d', 'once'))); % Remove subfolders that don't contain numbers
subfolder_names = {subfolders.name};
exp_dates = cellfun(@(x) regexp(x, '^\d{8}', 'match', 'once'), subfolder_names, 'UniformOutput', false);
exp_dates = exp_dates(~cellfun('isempty', exp_dates));

% Assign lambda and X_base based on exp_dates
lambda = zeros(1, length(exp_dates));
X_base = zeros(1, length(exp_dates));
flowfocusing_bandwidth = zeros(1, length(exp_dates));
for i = 1:length(exp_dates)
    switch exp_dates{i}
        case {'20231030', '20240415', '20241031'}
            lambda(i) = 30;
        case {'20240522', '20241111'}
            lambda(i) = 45;
        otherwise
            lambda(i) = NaN;
    end
    switch exp_dates{i}
        case {'20241031', '20241111'}
            X_base(i) = 6000;
            flowfocusing_bandwidth(i) = 100;
        case {'20231030', '20240415', '20240522'}
            X_base(i) = 9000;
            flowfocusing_bandwidth(i) = 150;
        otherwise
            X_base(i) = NaN;
    end
end

fiber_lengths = []; fiber_Ys = []; fiber_Xs = []; lambda_all = []; flowfocusing_bandwidth_all = [];
for exp_date_idx = [2,4,6,7] %1:length(subfolder_names)
    results_folder = fullfile(data_path, subfolder_names{exp_date_idx}, 'results');
    if isfolder(results_folder)
        mat_files = dir(fullfile(results_folder, '*.mat'));
        for mf = 1:length(mat_files)
            Y_base = str2double(extractBetween(mat_files(mf).name, 'Y', 'um'));
            load(fullfile(mat_files(mf).folder, mat_files(mf).name));
            for frame_idx = 1:improc
                for fiber_idx = 1 : FilNum
                    if ~exist('xy_after_selection', 'var') || isempty(xy_after_selection)
                        continue;
                    end
                    try 
                        fiber_L = xy_after_selection(fiber_idx).arclen_spl(frame_idx);
                        fiber_CoM = xy_after_selection(fiber_idx).centroid{frame_idx}*pixel_size;
                        % % % cartesian coordinate system (two flips in the following just to make the code robust)
                        fiber_CoM(2) = image_height - image_height + fiber_CoM(2);

                        fiber_lengths = [fiber_lengths; fiber_L*pixel_size];
                        fiber_Ys = [fiber_Ys; Y_base + fiber_CoM(2)];
                        fiber_Xs = [fiber_Xs; X_base(exp_date_idx) + fiber_CoM(1)];
                        lambda_all = [lambda_all; lambda(exp_date_idx)];
                        flowfocusing_bandwidth_all = [flowfocusing_bandwidth_all; flowfocusing_bandwidth(exp_date_idx)];
                        clear fiber_L fiber_CoM
                    catch
                        disp('No filament detected in this frame')
                    end
                end
            end
        end
    else
        exp2proc = dir(fullfile(data_path, subfolder_names{exp_date_idx}, 'Y*'));
        for ii = 1:length(exp2proc)
            Y_base = exp2proc(ii).name;
            Y_base = str2double(cell2mat((extractBetween(Y_base,'Y','um'))));
            Tracings = dir(fullfile(exp2proc(ii).folder, exp2proc(ii).name, 'Tracings_*.csv'));
            Vertices = dir(fullfile(exp2proc(ii).folder, exp2proc(ii).name, 'Vertices_*.csv'));

            for jj = 1:length(Tracings)
                Tracing_info = readmatrix(fullfile(Tracings(jj).folder, Tracings(jj).name),"NumHeaderLines",1);
                Vertice_info = readcell(fullfile(Vertices(jj).folder, Vertices(jj).name),"NumHeaderLines",1);

                fiber_lengths = [fiber_lengths; Tracing_info(:, 6)*pixel_size];
                fiber_index = Vertice_info(:,2);
                fiber_index = cell2mat(cellfun(@(x) str2double(x(2:end)), fiber_index, 'UniformOutput', false));
                for kk = 1:max(fiber_index)
                    fiber_xy = cell2mat(Vertice_info(fiber_index == kk,5:6));
                    fiber_CoM = mean(fiber_xy,1)*pixel_size;
                    % % % flip the fiber_y from image coordinate system to cartesian coordinate system
                    fiber_CoM(2) = image_height - fiber_CoM(2);

                    fiber_Ys = [fiber_Ys; Y_base + fiber_CoM(2)];
                    fiber_Xs = [fiber_Xs; X_base(exp_date_idx) + fiber_CoM(1)];
                    lambda_all = [lambda_all; lambda(exp_date_idx)];
                    flowfocusing_bandwidth_all = [flowfocusing_bandwidth_all; flowfocusing_bandwidth(exp_date_idx)];
                    clear fiber_xy fiber_CoM
                end
            end
        end
    end
end

% calculate the normalized fiber length and drift angle
fiber_lengths_normalized = fiber_lengths./lambda_all;
fiber_Ys_max = fiber_Ys;
fiber_Ys_min = fiber_Ys - flowfocusing_bandwidth_all;
fiber_Ys_min(fiber_Ys_min < 0) = 0; % set value to zeros if it's negative
beta_max_all = atand(fiber_Ys_max./fiber_Xs);
beta_min_all = atand(fiber_Ys_min./fiber_Xs);



%% Plotting: beta vs normalized fiber length (all data together)

% Separate data into two groups based on lambda
group1_idx = lambda_all == 30;
group2_idx = lambda_all == 45;

% Group 1: lambda = 30
fiber_lengths_g1 = fiber_lengths(group1_idx);
fiber_lengths_normalized_g1 = fiber_lengths_normalized(group1_idx);
fiber_Ys_max_g1 = fiber_Ys_max(group1_idx);
fiber_Xs_g1 = fiber_Xs(group1_idx);
beta_max_all_g1 = beta_max_all(group1_idx);

% Group 2: lambda = 45
fiber_lengths_g2 = fiber_lengths(group2_idx);
fiber_lengths_normalized_g2 = fiber_lengths_normalized(group2_idx);
fiber_Ys_max_g2 = fiber_Ys_max(group2_idx);
fiber_Xs_g2 = fiber_Xs(group2_idx);
beta_max_all_g2 = beta_max_all(group2_idx);

% Bin and average for group 1 (lambda = 30)
bin_size = 0.12; bin_number = 20;
L_mean_g1 = zeros(1, bin_number+1);
L_std_g1 = zeros(1, bin_number+1);
Y_max_mean_g1 = zeros(1, bin_number+1);
Y_max_std_g1 = zeros(1, bin_number+1);
X_mean_g1 = zeros(1, bin_number+1);
for ii = 0:bin_number
    if ii < bin_number
        idx = fiber_lengths_normalized_g1 > ii*bin_size & fiber_lengths_normalized_g1 < (ii+1)*bin_size;
    else
        idx = fiber_lengths_normalized_g1 > ii*bin_size;
    end
    L_mean_g1(ii+1) = mean(fiber_lengths_normalized_g1(idx));
    L_std_g1(ii+1) = std(fiber_lengths_normalized_g1(idx));
    Y_max_mean_g1(ii+1) = mean(fiber_Ys_max_g1(idx));
    Y_max_std_g1(ii+1) = std(fiber_Ys_max_g1(idx));
    X_mean_g1(ii+1) = mean(fiber_Xs_g1(idx));
end

% Bin and average for group 2 (lambda = 45)
L_mean_g2 = zeros(1, bin_number+1);
L_std_g2 = zeros(1, bin_number+1);
Y_max_mean_g2 = zeros(1, bin_number+1);
Y_max_std_g2 = zeros(1, bin_number+1);
X_mean_g2 = zeros(1, bin_number+1);
for ii = 0:bin_number
    if ii < bin_number
        idx = fiber_lengths_normalized_g2 > ii*bin_size & fiber_lengths_normalized_g2 < (ii+1)*bin_size;
    else
        idx = fiber_lengths_normalized_g2 > ii*bin_size;
    end
    L_mean_g2(ii+1) = mean(fiber_lengths_normalized_g2(idx));
    L_std_g2(ii+1) = std(fiber_lengths_normalized_g2(idx));
    Y_max_mean_g2(ii+1) = mean(fiber_Ys_max_g2(idx));
    Y_max_std_g2(ii+1) = std(fiber_Ys_max_g2(idx));
    X_mean_g2(ii+1) = mean(fiber_Xs_g2(idx));
end

% Plotting:
figure('Color','w','Position',[100 100 800 600]);

% plot the original data points: group 1 (lambda = 30)
plot(fiber_lengths_normalized_g1, beta_max_all_g1, 'color', [203, 213, 232]/255, ...
    'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 8, 'LineWidth', 1.2);
hold on;
% plot the original data points: group 2 (lambda = 45)
plot(fiber_lengths_normalized_g2, beta_max_all_g2, 'color', [253, 205, 172]/255, ...
    'LineStyle', 'none', 'Marker', '^', 'MarkerSize', 8, 'LineWidth', 1.2);
hold on;

% errorbar plot: bin average based on normalized fiber length
% errorbar plot: group 1 (lambda = 30)
errorbar(L_mean_g1, atand(Y_max_mean_g1./X_mean_g1), atand(Y_max_std_g1./X_mean_g1), ...
    atand(Y_max_std_g1./X_mean_g1), L_std_g1, L_std_g1, 'color', [117, 112, 179]/255, ...
    'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
% errorbar plot: group 2 (lambda = 45)
errorbar(L_mean_g2, atand(Y_max_mean_g2./X_mean_g2), atand(Y_max_std_g2./X_mean_g2), ...
    atand(Y_max_std_g2./X_mean_g2), L_std_g2, L_std_g2, 'color', [217, 95, 2]/255, ...
    'LineStyle', 'none', 'Marker', '^', 'MarkerSize', 10, 'LineWidth', 2);
hold on;

ylabel('$\beta\ (\mathrm{^\circ})$','FontSize', 24,'Interpreter', 'latex');
xlabel('$L/\lambda$','FontSize', 24,'Interpreter', 'latex');
xlim([0 3.3]); ylim([0 10])
set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 22)

legend({'', '', '$\lambda=30\mathrm{\mu m}$', ...
    '$\lambda=45\mathrm{\mu m}$'},'Interpreter', 'latex', 'FontSize', 20)

% text(2.6, 0.7, '$\lambda/R=3$', 'FontSize', 24, 'FontWeight', 'bold', ...
%     'Interpreter','latex','BackgroundColor','w');

f=gcf;
saveas(f, ['Z:\Processing & Results\Actin Filaments in Porous Media\Actin-DLD\' ...
    'Figures\Beta_Length_errorbarPlot_ALLexp_BinOnNormalized.fig'])

f=gcf;
exportgraphics(f,['Z:\Processing & Results\Actin Filaments in Porous Media\Actin-DLD\' ...
    'Figures\Beta_Length_errorbarPlot_ALLexp_BinOnNormalized.png'],'Resolution',100)

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['Z:\Processing & Results\Actin Filaments in Porous Media\Actin-DLD\' ...
    'Figures\Beta_Length_errorbarPlot_ALLexp_BinOnNormalized.eps']);



%% Plotting: sorting efficiency

% Separate data into two groups based on X_base
% Here use flowfocusing_bandwidth to separate the two groups since X_base and flowfocusing_bandwidth are coupled in this dataset
group_Xbase_6000_idx = flowfocusing_bandwidth_all == 100;
group_Xbase_9000_idx = flowfocusing_bandwidth_all == 150;

% Group 1: X_base = 6000
fiber_lengths_normalized_Xbase6000 = fiber_lengths_normalized(group_Xbase_6000_idx);
fiber_Ys_max_Xbase6000 = fiber_Ys_max(group_Xbase_6000_idx);
lambda_all_Xbase6000 = lambda_all(group_Xbase_6000_idx);

% Group 2: X_base = 9000
fiber_lengths_normalized_Xbase9000 = fiber_lengths_normalized(group_Xbase_9000_idx);
fiber_Ys_max_Xbase9000 = fiber_Ys_max(group_Xbase_9000_idx);
lambda_all_Xbase9000 = lambda_all(group_Xbase_9000_idx);

% Loop to plot for both X_base = 6000 and X_base = 9000

% (1) Method 1: Histogram counts along Y for different fiber lengths at given X_base
for Xbase_val = [6000, 9000]
    if Xbase_val == 6000
        fiber_lengths_normalized_sorting_efficiency_plot = fiber_lengths_normalized_Xbase6000;
        fiber_Y_normalized_sorting_efficiency_plot = fiber_Ys_max_Xbase6000 ./ lambda_all_Xbase6000;
        suffix = '_Xbase6000';
    else
        fiber_lengths_normalized_sorting_efficiency_plot = fiber_lengths_normalized_Xbase9000;
        fiber_Y_normalized_sorting_efficiency_plot = fiber_Ys_max_Xbase9000 ./ lambda_all_Xbase9000;
        suffix = '_Xbase9000';
    end


    % Group the data based on normalized fiber length
    short_idx = fiber_lengths_normalized_sorting_efficiency_plot < 1.3;
    middle_idx = fiber_lengths_normalized_sorting_efficiency_plot >= 1.3 & fiber_lengths_normalized_sorting_efficiency_plot <= 1.8;
    long_idx = fiber_lengths_normalized_sorting_efficiency_plot > 1.8;

    short_fiber_Y = fiber_Y_normalized_sorting_efficiency_plot(short_idx);
    short_fiber_length = fiber_lengths_normalized_sorting_efficiency_plot(short_idx);

    middle_fiber_Y = fiber_Y_normalized_sorting_efficiency_plot(middle_idx);
    middle_fiber_length = fiber_lengths_normalized_sorting_efficiency_plot(middle_idx);

    long_fiber_Y = fiber_Y_normalized_sorting_efficiency_plot(long_idx);
    long_fiber_length = fiber_lengths_normalized_sorting_efficiency_plot(long_idx);

    % Plot histograms for each fiber length category
    figure('Color','w','Position',[0 100 800 600]);
    hold on;
    histogram(short_fiber_Y, 'BinWidth', 2, 'FaceAlpha', 0.6, 'FaceColor', [255, 18, 18]/255, 'EdgeColor', 'none');
    histogram(middle_fiber_Y, 'BinWidth', 2, 'FaceAlpha', 0.6, 'FaceColor', [12, 178, 34]/255, 'EdgeColor', 'none');
    histogram(long_fiber_Y, 'BinWidth', 2, 'FaceAlpha', 0.6, 'FaceColor', [14, 35, 202]/255, 'EdgeColor', 'none');
    xlabel('$y/\lambda$','FontSize', 24,'Interpreter', 'latex');
    ylabel('Fiber count','FontSize', 24,'Interpreter', 'latex');
    legend({'Short ($L/\lambda < 1.3$)', 'Medium ($1.3 \leq L/\lambda \leq 1.8$)', 'Long ($L/\lambda > 1.8$)'}, 'Interpreter', 'latex', 'FontSize', 20);
    set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24)
    hold off;

    f=gcf;
    exportgraphics(f,['Z:\Processing & Results\Actin Filaments in Porous Media\Actin-DLD\' ...
        'Figures\Sorting_efficiency_counts_along_Y' suffix '.png'],'Resolution',100)

    f=gcf;
    exportgraphics(f,['Z:\Processing & Results\Actin Filaments in Porous Media\Actin-DLD\' ...
        'Figures\Sorting_efficiency_counts_along_Y' suffix '.pdf'],'ContentType','vector')

    % Normalize histograms by total count for each fiber type
    figure('Color','w','Position',[900 100 800 600]);
    hold on;
    histogram(short_fiber_Y, 'BinWidth', 2, 'FaceAlpha', 0.6, 'FaceColor', [255, 18, 18]/255, ...
        'EdgeColor', 'none', 'Normalization', 'probability');
    histogram(middle_fiber_Y, 'BinWidth', 2, 'FaceAlpha', 0.6, 'FaceColor', [12, 178, 34]/255, ...
        'EdgeColor', 'none', 'Normalization', 'probability');
    histogram(long_fiber_Y, 'BinWidth', 2, 'FaceAlpha', 0.6, 'FaceColor', [14, 35, 202]/255, ...
        'EdgeColor', 'none', 'Normalization', 'probability');

    xlabel('$y/\lambda$','FontSize', 24,'Interpreter', 'latex');
    ylabel('Fraction of fibers','FontSize', 24,'Interpreter', 'latex');
    legend({'Short', 'Medium', 'Long'}, 'Interpreter', 'latex', 'FontSize', 20);
    set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24)
    hold off;

    f=gcf;
    exportgraphics(f,['Z:\Processing & Results\Actin Filaments in Porous Media\Actin-DLD\' ...
        'Figures\Sorting_efficiency_counts_along_Y_normalized' suffix '.png'],'Resolution',100)

    f=gcf;
    exportgraphics(f,['Z:\Processing & Results\Actin Filaments in Porous Media\Actin-DLD\' ...
        'Figures\Sorting_efficiency_counts_along_Y_normalized' suffix '.pdf'],'ContentType','vector')

    % Plot smoothed curves instead of histograms for each fiber length category
    figure('Color','w','Position',[1800 100 800 600]);
    hold on;

    % Define bin edges for normalized Y position
    edges = 0:2:max(fiber_Y_normalized_sorting_efficiency_plot)+1;

    % Short fibers
    [counts_short, ~] = histcounts(short_fiber_Y, edges, 'Normalization', 'probability');
    centers_short = edges(1:end-1) + diff(edges)/2;
    plot(centers_short, smooth(counts_short, 3), '-','LineWidth',2,'Color',[255, 18, 18]/255);

    % Medium fibers
    [counts_middle, ~] = histcounts(middle_fiber_Y, edges, 'Normalization', 'probability');
    centers_middle = edges(1:end-1) + diff(edges)/2;
    plot(centers_middle, smooth(counts_middle, 3), '-','LineWidth',2,'Color',[12, 178, 34]/255);

    % Long fibers
    [counts_long, ~] = histcounts(long_fiber_Y, edges, 'Normalization', 'probability');
    centers_long = edges(1:end-1) + diff(edges)/2;
    plot(centers_long, smooth(counts_long, 3), '-','LineWidth',2,'Color',[14, 35, 202]/255);

    xlabel('$y/\lambda$','FontSize', 24,'Interpreter', 'latex');
    ylabel('Fraction of fibers','FontSize', 24,'Interpreter', 'latex');
    legend({'Short', 'Medium', 'Long'}, 'Interpreter', 'latex', 'FontSize', 20);
    set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24)
    hold off;

    f=gcf;
    exportgraphics(f,['Z:\Processing & Results\Actin Filaments in Porous Media\Actin-DLD\' ...
        'Figures\Sorting_efficiency_counts_along_Y_normalized_curve' suffix '.png'],'Resolution',100)
end


% (2) Method 2: Histogram counts of normalized fiber length at different Y sections (two sections)
for Xbase_val = [6000, 9000]
    if Xbase_val == 6000
        fiber_lengths_normalized_sorting_efficiency_plot = fiber_lengths_normalized_Xbase6000;
        fiber_Y_normalized_sorting_efficiency_plot = fiber_Ys_max_Xbase6000 ./ lambda_all_Xbase6000;
        suffix = '_Xbase6000';
        section_bounds = 20;
    else
        fiber_lengths_normalized_sorting_efficiency_plot = fiber_lengths_normalized_Xbase9000;
        fiber_Y_normalized_sorting_efficiency_plot = fiber_Ys_max_Xbase9000 ./ lambda_all_Xbase9000;
        suffix = '_Xbase9000';
        section_bounds = 20;
    end

    % Define three Y sections using the specified boundaries
    Y_sections_1_idx = fiber_Y_normalized_sorting_efficiency_plot <= section_bounds(1);
    Y_sections_2_idx = fiber_Y_normalized_sorting_efficiency_plot > section_bounds(1);

    section_1_lengths = fiber_lengths_normalized_sorting_efficiency_plot(Y_sections_1_idx);
    section_2_lengths = fiber_lengths_normalized_sorting_efficiency_plot(Y_sections_2_idx);        

    % Determine common x-axis limits based on all sections
    all_lengths = [section_1_lengths; section_2_lengths];
    x_min = min(all_lengths);
    x_max = max(all_lengths);

    figure('Color','w','Position',[100 100 700 600]);
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    % Section 2
    nexttile;
    histogram(section_2_lengths, 'BinWidth', 0.2, 'FaceAlpha', 0.7, 'FaceColor', [0, 158, 115]/255, 'EdgeColor', 'none');
    ylabel('Fiber count', 'FontSize', 24, 'Interpreter', 'latex');
    title(sprintf('$y/\\lambda > %.0f$', section_bounds(1)), 'FontSize', 18, 'Interpreter', 'latex');
    set(gca, 'Box', 'On', 'FontSize', 20, 'XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5);
    set(gca, 'XTickLabel', []);
    xlim([x_min x_max]);

    % Section 1
    nexttile;
    histogram(section_1_lengths, 'BinWidth', 0.2, 'FaceAlpha', 0.7, 'FaceColor', [230, 159, 0]/255, 'EdgeColor', 'none');
    xlabel('$L/\lambda$', 'FontSize', 20, 'Interpreter', 'latex');
    ylabel('Fiber count', 'FontSize', 24, 'Interpreter', 'latex');
    title(sprintf('$y/\\lambda \\leq %.0f$', section_bounds(1)), 'FontSize', 18, 'Interpreter', 'latex');
    set(gca, 'Box', 'On', 'FontSize', 20, 'XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5);
    xlim([x_min x_max]);

    f=gcf;
    exportgraphics(f,['Z:\Processing & Results\Actin Filaments in Porous Media\Actin-DLD\' ...
        'Figures\Sorting_efficiency_counts_at_diff_Y_two_sections' suffix '.png'],'Resolution',100)

    f=gcf;
    exportgraphics(f,['Z:\Processing & Results\Actin Filaments in Porous Media\Actin-DLD\' ...
        'Figures\Sorting_efficiency_counts_at_diff_Y_two_sections' suffix '.pdf'],'ContentType','vector');
end


% (2) Method 2: Histogram counts of normalized fiber length at different Y sections
for Xbase_val = [6000, 9000]
    if Xbase_val == 6000
        fiber_lengths_normalized_sorting_efficiency_plot = fiber_lengths_normalized_Xbase6000;
        fiber_Y_normalized_sorting_efficiency_plot = fiber_Ys_max_Xbase6000 ./ lambda_all_Xbase6000;
        suffix = '_Xbase6000';
        section_bounds = [15, 30];
    else
        fiber_lengths_normalized_sorting_efficiency_plot = fiber_lengths_normalized_Xbase9000;
        fiber_Y_normalized_sorting_efficiency_plot = fiber_Ys_max_Xbase9000 ./ lambda_all_Xbase9000;
        suffix = '_Xbase9000';
        section_bounds = [18, 36];
    end

    % Define three Y sections using the specified boundaries
    Y_sections_1_idx = fiber_Y_normalized_sorting_efficiency_plot <= section_bounds(1);
    Y_sections_2_idx = fiber_Y_normalized_sorting_efficiency_plot > section_bounds(1) & fiber_Y_normalized_sorting_efficiency_plot <= section_bounds(2);
    Y_sections_3_idx = fiber_Y_normalized_sorting_efficiency_plot > section_bounds(2);

    section_1_lengths = fiber_lengths_normalized_sorting_efficiency_plot(Y_sections_1_idx);
    section_2_lengths = fiber_lengths_normalized_sorting_efficiency_plot(Y_sections_2_idx);        
    section_3_lengths = fiber_lengths_normalized_sorting_efficiency_plot(Y_sections_3_idx);

    % Determine common x-axis limits based on all sections
    all_lengths = [section_1_lengths; section_2_lengths; section_3_lengths];
    x_min = min(all_lengths);
    x_max = max(all_lengths);

    figure('Color','w','Position',[100 100 700 900]);
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    % Section 3
    nexttile;
    histogram(section_3_lengths, 'BinWidth', 0.2, 'FaceAlpha', 0.7, 'FaceColor', [0, 158, 115]/255, 'EdgeColor', 'none');
    ylabel('Fiber count', 'FontSize', 24, 'Interpreter', 'latex');
    title(sprintf('$y/\\lambda > %.0f$', section_bounds(2)), 'FontSize', 18, 'Interpreter', 'latex');
    set(gca, 'Box', 'On', 'FontSize', 20, 'XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5);
    set(gca, 'XTickLabel', []);
    xlim([x_min x_max]);

    % Section 2
    nexttile;
    histogram(section_2_lengths, 'BinWidth', 0.2, 'FaceAlpha', 0.7, 'FaceColor', [230, 159, 0]/255, 'EdgeColor', 'none');
    ylabel('Fiber count', 'FontSize', 24, 'Interpreter', 'latex');
    title(sprintf('$%.0f < y/\\lambda \\leq %.0f$', section_bounds(1), section_bounds(2)), 'FontSize', 18, 'Interpreter', 'latex');
    set(gca, 'Box', 'On', 'FontSize', 20, 'XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5);
    set(gca, 'XTickLabel', []);
    xlim([x_min x_max]);

    % Section 1
    nexttile;
    histogram(section_1_lengths, 'BinWidth', 0.2, 'FaceAlpha', 0.7, 'FaceColor', [86, 180, 233]/255, 'EdgeColor', 'none');
    xlabel('$L/\lambda$', 'FontSize', 20, 'Interpreter', 'latex');
    ylabel('Fiber count', 'FontSize', 24, 'Interpreter', 'latex');
    title(sprintf('$y/\\lambda \\leq %.0f$', section_bounds(1)), 'FontSize', 18, 'Interpreter', 'latex');
    set(gca, 'Box', 'On', 'FontSize', 20, 'XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5);
    xlim([x_min x_max]);

    f=gcf;
    exportgraphics(f,['Z:\Processing & Results\Actin Filaments in Porous Media\Actin-DLD\' ...
        'Figures\Sorting_efficiency_counts_at_diff_Y' suffix '.png'],'Resolution',100)

    f=gcf;
    exportgraphics(f,['Z:\Processing & Results\Actin Filaments in Porous Media\Actin-DLD\' ...
        'Figures\Sorting_efficiency_counts_at_diff_Y' suffix '.pdf'],'ContentType','vector');
end


% (2) Method 2: Histogram counts of normalized fiber length at different Y sections (four sections)
for Xbase_val = [6000, 9000]
    if Xbase_val == 6000
        fiber_lengths_normalized_sorting_efficiency_plot = fiber_lengths_normalized_Xbase6000;
        fiber_Y_normalized_sorting_efficiency_plot = fiber_Ys_max_Xbase6000 ./ lambda_all_Xbase6000;
        suffix = '_Xbase6000_four_sections';
        section_bounds = [10, 20, 30];
    else
        fiber_lengths_normalized_sorting_efficiency_plot = fiber_lengths_normalized_Xbase9000;
        fiber_Y_normalized_sorting_efficiency_plot = fiber_Ys_max_Xbase9000 ./ lambda_all_Xbase9000;
        suffix = '_Xbase9000_four_sections';
        section_bounds = [9, 24, 39];
    end

    % Define four Y sections using the specified boundaries
    Y_sections_1_idx = fiber_Y_normalized_sorting_efficiency_plot <= section_bounds(1);
    Y_sections_2_idx = fiber_Y_normalized_sorting_efficiency_plot > section_bounds(1) & fiber_Y_normalized_sorting_efficiency_plot <= section_bounds(2);
    Y_sections_3_idx = fiber_Y_normalized_sorting_efficiency_plot > section_bounds(2) & fiber_Y_normalized_sorting_efficiency_plot <= section_bounds(3);
    Y_sections_4_idx = fiber_Y_normalized_sorting_efficiency_plot > section_bounds(3);

    section_1_lengths = fiber_lengths_normalized_sorting_efficiency_plot(Y_sections_1_idx);
    section_2_lengths = fiber_lengths_normalized_sorting_efficiency_plot(Y_sections_2_idx);
    section_3_lengths = fiber_lengths_normalized_sorting_efficiency_plot(Y_sections_3_idx);
    section_4_lengths = fiber_lengths_normalized_sorting_efficiency_plot(Y_sections_4_idx);

    % Determine common x-axis limits based on all sections
    all_lengths = [section_1_lengths; section_2_lengths; section_3_lengths; section_4_lengths];
    x_min = min(all_lengths);
    x_max = max(all_lengths);

    figure('Color','w','Position',[100 100 700 1200]);
    tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

    % Section 4
    nexttile;
    histogram(section_4_lengths, 'BinWidth', 0.2, 'FaceAlpha', 0.7, 'FaceColor', [0, 158, 115]/255, 'EdgeColor', 'none');
    ylabel('Fiber count', 'FontSize', 16, 'Interpreter', 'latex');
    title(sprintf('$y/\\lambda > %.0f$', section_bounds(3)), 'FontSize', 18, 'Interpreter', 'latex');
    set(gca, 'Box', 'On', 'FontSize', 14, 'XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5);
    set(gca, 'XTickLabel', []);
    xlim([x_min x_max]);

    % Section 3
    nexttile;
    histogram(section_3_lengths, 'BinWidth', 0.2, 'FaceAlpha', 0.7, 'FaceColor', [230, 159, 0]/255, 'EdgeColor', 'none');
    ylabel('Fiber count', 'FontSize', 16, 'Interpreter', 'latex');
    title(sprintf('$%.0f < y/\\lambda \\leq %.0f$', section_bounds(2), section_bounds(3)), 'FontSize', 18, 'Interpreter', 'latex');
    set(gca, 'Box', 'On', 'FontSize', 14, 'XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5);
    set(gca, 'XTickLabel', []);
    xlim([x_min x_max]);

    % Section 2
    nexttile;
    histogram(section_2_lengths, 'BinWidth', 0.2, 'FaceAlpha', 0.7, 'FaceColor', [255, 140, 0]/255, 'EdgeColor', 'none');
    ylabel('Fiber count', 'FontSize', 16, 'Interpreter', 'latex');
    title(sprintf('$%.0f < y/\\lambda \\leq %.0f$', section_bounds(1), section_bounds(2)), 'FontSize', 18, 'Interpreter', 'latex');
    set(gca, 'Box', 'On', 'FontSize', 14, 'XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5);
    set(gca, 'XTickLabel', []);
    xlim([x_min x_max]);

    % Section 1
    nexttile;
    histogram(section_1_lengths, 'BinWidth', 0.2, 'FaceAlpha', 0.7, 'FaceColor', [86, 180, 233]/255, 'EdgeColor', 'none');
    xlabel('$L/\lambda$', 'FontSize', 20, 'Interpreter', 'latex');
    ylabel('Fiber count', 'FontSize', 16, 'Interpreter', 'latex');
    title(sprintf('$y/\\lambda \\leq %.0f$', section_bounds(1)), 'FontSize', 18, 'Interpreter', 'latex');
    set(gca, 'Box', 'On', 'FontSize', 14, 'XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5);
    xlim([x_min x_max]);

    f=gcf;
    exportgraphics(f,['Z:\Processing & Results\Actin Filaments in Porous Media\Actin-DLD\' ...
        'Figures\Sorting_efficiency_counts_at_diff_Y_four_sections' suffix '.png'],'Resolution',100)

    f=gcf;
    exportgraphics(f,['Z:\Processing & Results\Actin Filaments in Porous Media\Actin-DLD\' ...
        'Figures\Sorting_efficiency_counts_at_diff_Y_four_sections' suffix '.pdf'],'ContentType','vector');
end