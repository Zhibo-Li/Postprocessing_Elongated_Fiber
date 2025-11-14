clear; close all; clc;

% Prompt the user to input parameters
pixel_size = 0.1; % the magnification of the objective (um/pixel)
Q = 10; % the flow rate in nL/s

image_height = 2048; % the height of the image, unit: pixel

% Precompute the colormap for whole fiber colorcoded plot
TheColor_whole_fiber_1 = cmocean('-ice', 18);
TheColor_whole_fiber_2 = cmocean('solar', 18);
TheColor_whole_fiber = [repmat([0.9 0.9 0.9], 75, 1); TheColor_whole_fiber_1(end-15:end-1, :); ...
    TheColor_whole_fiber_2(2:16, :); repmat([0.9 0.9 0.9], 75, 1)];

% Set the path for loading data
exp2proc = 'Z:\Processing & Results\ActinSuspUnderFlow\FiberAngleUnderShear_W150H500\';

% Choose one or multiple files to load
Selected_Names = uigetfile_n_dir(exp2proc, 'Choose the data to be processed:');
% Get new app 'uigetfile_n_dir', and the follows are needed:
% 1. Download MATLAB compiler SDK.
% 2. Rename the function to "uigetfile_n_dir" in line 1.
% 3. Change line 6: "start_path == ''" to "isempty(start_path)".
% 4. Remove "|| start_path == 0" in line 6.


ExpDates = length(Selected_Names);
% Initialize Fiber_info struct array (one element per selected experiment)
Fiber_info = repmat(struct(...
    'fiber_spl', [], ...                  % cell array of fiber splines for this experiment
    'susp_conc', [], ...                  % suspension concentration (%)
    'Z', [], ...                          % Z position (um)
    'L0_all', [], ...                     % fiber lengths (um) for this experiment
    'Chi_all', [], ...                    % fiber orientations (deg) for this experiment
    'CoM_y_all_normalized', [], ...       % normalized CoM y-positions for this experiment
    'S', [], ...                          % order parameter for this experiment
    'Q', Q), 1, 100);
n = 1;
for ii = 1:ExpDates

    % read the txt file
    TheBackground_boundary = readmatrix(fullfile(Selected_Names{ii}, 'TheBackground.txt')); % preferred for numeric text
    Lower_boundary = TheBackground_boundary(2); % in pixel
    Upper_boundary = TheBackground_boundary(1); % in pixel

    matPattern = fullfile(Selected_Names{ii}, 'AfterAveBGR', 'results', '*.mat');
    matFiles = dir(matPattern);
    if isempty(matFiles)
        warning('No .mat files found: %s', matPattern);
    end
    for mf = 1:numel(matFiles)
        matFilePath = fullfile(matFiles(mf).folder, matFiles(mf).name);
        
        [~, theFILE] = fileparts(matFilePath);
        Z = str2double(extractBetween(theFILE,'Z','um'));
        concStr = extractBetween(theFILE, 'X-Mid_','Percent');
        concStr = replace(concStr, 'o', '.');
        susp_conc = str2double(concStr);

        fprintf('Loading: %s\n', matFilePath);
        load(matFilePath, 'xy_after_selection', 'FilNum', 'improc'); % loads variables from this .mat into the workspace
        L0_all = []; % fiber lengths (um)
        Chi_all = []; % fiber orientations (deg)
        CoM_y_all = []; % fiber CoM y-positions (pixel)
        fiber_spl = cell(0);
        for frame_idx = 1:improc
            for fiber_idx = 1 : FilNum
                try
                    spl = xy_after_selection(fiber_idx).spl{frame_idx};
                    CoM_xy = xy_after_selection(fiber_idx).centroid{frame_idx};
                    % % % cartesian coordinate system (two flips in the following just to make the code robust)
                    CoM_xy(2) = image_height - image_height + CoM_xy(2);

                    % Orientation of whole fiber
                    Gyr = 1/size(spl,1) * [sum((spl(:, 1)-CoM_xy(1)).^2),  sum((spl(:, 1)-CoM_xy(1)) .* (spl(:, 2)-CoM_xy(2)));
                        sum((spl(:, 2)-CoM_xy(2)) .* (spl(:, 1)-CoM_xy(1))), sum((spl(:, 2)-CoM_xy(2)).^2)];

                    [eigenV,eigenD] = eig(Gyr);
                    [d,ind] = sort(diag(eigenD));
                    Ds = eigenD(ind,ind);
                    Vs = eigenV(:,ind);
                    Chi = atand(Vs(2,2)/Vs(1,2));

                    % % %% Create a colorcoded plot (color based on the orientation of the whole fiber)
                    % % color_indices_wholefiber = floor(Chi+90) + 1;
                    % % plot(spl(:,1), spl(:,2), 'Color', TheColor_whole_fiber(color_indices_wholefiber, :), 'LineWidth', 3);
                    % % hold on
                    % % xlim([0 2024]); ylim([0 2024])

                    spl_Ls = xy_after_selection(fiber_idx).arclen_spl(frame_idx);
                    fiber_spl = [fiber_spl; {spl}];
                    L0_all = [L0_all; spl_Ls * pixel_size];
                    Chi_all = [Chi_all; Chi];
                    CoM_y_all = [CoM_y_all; CoM_xy(2)];
                catch
                    disp('No filament detected in this frame')
                end
            end
        end
        
        % Order parameter
        S = 0.5*(3*mean(cosd(Chi_all).^2)-1);

        Fiber_info(n).fiber_spl = fiber_spl;
        Fiber_info(n).susp_conc = susp_conc;
        Fiber_info(n).Z = Z;
        Fiber_info(n).L0_all = L0_all;
        Fiber_info(n).Chi_all = Chi_all;
        CoM_y_all_normalized = (CoM_y_all - Lower_boundary)/(Upper_boundary - Lower_boundary);
        Fiber_info(n).CoM_y_all_normalized = CoM_y_all_normalized;
        Fiber_info(n).S = S;
        % Fiber_info(n).Q = Q;
        clear L0_all Chi_all CoM_y_all
        
        n = n + 1;
    end
end

% remove the empty elements in Fiber_info
Fiber_info = Fiber_info(~cellfun('isempty', {Fiber_info.Z}));

% Group data by suspension concentration and Z position
[unique_pairs, ~, group_idx] = unique([[Fiber_info.susp_conc]', [Fiber_info.Z]'], 'rows');
num_groups = size(unique_pairs, 1);

% Combine data for each group (same susp_conc and Z)
Combined_info = repmat(struct(...
    'fiber_spl', [], ...
    'susp_conc', [], ...
    'Z', [], ...
    'L0_all', [], ...
    'Chi_all', [], ...
    'CoM_y_all_normalized', [], ...
    'S', [], ...
    'Q', []), 1, num_groups);

for g = 1:num_groups
    % Find all entries belonging to this group
    group_indices = find(group_idx == g);
    
    % Combine all data from this group
    combined_fiber_spl = {};
    combined_L0 = [];
    combined_Chi = [];
    combined_CoM_y = [];
    
    for idx = group_indices'
        combined_fiber_spl = [combined_fiber_spl; Fiber_info(idx).fiber_spl];
        combined_L0 = [combined_L0; Fiber_info(idx).L0_all];
        combined_Chi = [combined_Chi; Fiber_info(idx).Chi_all];
        combined_CoM_y = [combined_CoM_y; Fiber_info(idx).CoM_y_all_normalized];
    end
    
    % Recalculate order parameter for combined data
    combined_S = 0.5*(3*mean(cosd(combined_Chi).^2)-1);
    
    % Store combined data
    Combined_info(g).fiber_spl = combined_fiber_spl;
    Combined_info(g).susp_conc = unique_pairs(g, 1);
    Combined_info(g).Z = unique_pairs(g, 2);
    Combined_info(g).L0_all = combined_L0;
    Combined_info(g).Chi_all = combined_Chi;
    Combined_info(g).CoM_y_all_normalized = combined_CoM_y;
    Combined_info(g).S = combined_S;
    Combined_info(g).Q = Q;
end


%% Create a colorcoded plot (color based on the orientation of the whole fiber)
% Create a colorcoded plot for each combined condition
for gg = 1:num_groups
    
    % Get data for this group
    fiber_spl = Combined_info(gg).fiber_spl;
    Chi_all = Combined_info(gg).Chi_all;
    susp_conc = Combined_info(gg).susp_conc;
    Z_pos = Combined_info(gg).Z;
    CoM_y_all_normalized = Combined_info(gg).CoM_y_all_normalized;
    
    % plot the colorcoded fibers
    figure("Color", "white", "Position", [100, 100, 900, 900]);

    % Plot all fibers in this group
    for f = 1:length(fiber_spl)
        spl = fiber_spl{f};
        Chi = Chi_all(f);
        color_indices_wholefiber = floor(Chi+90) + 1;
        plot(spl(:,1), spl(:,2), 'Color', TheColor_whole_fiber(color_indices_wholefiber, :), 'LineWidth', 2);
        hold on
    end
    xlim([0 2024]); ylim([0 2024])
    % Add colorbar for TheColor_whole_fiber
    colormap(TheColor_whole_fiber);
    clim([1 size(TheColor_whole_fiber, 1)]);
    set(gcf, "Color", "white", "Position", [100, 100, 700, 600]);
    c = colorbar('Ticks', linspace(1, size(TheColor_whole_fiber, 1), 5), ...
                 'TickLabels', round(linspace(-90, 90, 5)), ...
                 'FontSize', 14, 'TickLabelInterpreter', 'LaTeX');
    c.Label.String = 'Orientation $\theta$ (deg)';
    c.Label.Interpreter = 'LaTeX';
    title(sprintf('Concentration: %.2f%%, Z: %.1f um, S: %.3f', susp_conc, Z_pos, Combined_info(gg).S));
    xlabel('X');
    ylabel('Y');
    axis equal; axis off
    
    % Save the figure
    saveas(gcf, fullfile(exp2proc, sprintf('ColorCodedFibers_Conc%.2f_Z%.1f.png', susp_conc, Z_pos)));


    

    
    % Plot PDF of Chi_all and fit with Gaussian
    figure("Color", "white", "Position", [100, 100, 800, 600]);
    histogram(Chi_all, 'Normalization', 'pdf', 'BinWidth', 1, 'FaceAlpha', 0.6, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
    hold on;

    % Robust Gaussian fit ignoring heavy tail (MAD-based trimming)
    Chi = Chi_all(~isnan(Chi_all));
    if isempty(Chi)
        pd.mu = NaN; pd.sigma = NaN;
        x_fit = []; y_fit = [];
    else
        if numel(Chi) >= 5
            mu0 = median(Chi);
            s0  = 1.4826*mad(Chi, 1);              % robust scale ~ std
            if s0 <= eps, s0 = max(std(Chi,1), eps); end
            k = 2.5;                                % reject tails beyond ~2.5 MADs
            inliers = abs(Chi - mu0) <= k*s0;
            if nnz(inliers) < max(5, round(0.5*numel(Chi)))
                inliers = abs(Chi - mu0) <= 3.5*s0; % relax if too few inliers
            end
            Chi_fit = Chi(inliers);
        else
            Chi_fit = Chi;
        end
        pd.mu = mean(Chi_fit);
        pd.sigma = std(Chi_fit, 1);
        if ~isfinite(pd.sigma) || pd.sigma <= 0, pd.sigma = max(std(Chi,1), 1e-6); end
        x_fit = linspace(min(Chi), max(Chi), 400);
        y_fit = (1./(pd.sigma*sqrt(2*pi))) .* exp(-0.5*((x_fit - pd.mu)./pd.sigma).^2);
    end
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

    xlabel('Orientation $\theta$ (deg)', 'Interpreter', 'latex');
    ylabel('Probability Density', 'Interpreter', 'latex');
    title(sprintf(['Actin concentration: %.2f\\%%, Z: %.1f um\n' ...
                   '$\\mu = %.2f^{\\circ},\\ \\sigma = %.2f^{\\circ}$'], ...
                 susp_conc, Z_pos, pd.mu, pd.sigma), 'Interpreter', 'latex');
    legend('Data', sprintf('Gaussian Fit'), 'Location', 'best', 'Interpreter', 'latex');
    grid on;
    xlim([-90 90]);
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');

    % Save the figure
    saveas(gcf, fullfile(exp2proc, sprintf('Chi_PDF_Fit_Conc%.2f_Z%.1f.png', susp_conc, Z_pos)));





    % Split by CoM_y_all_normalized threshold and plot both groups together
    idx_low  = CoM_y_all_normalized < 0.5;
    idx_high = CoM_y_all_normalized >= 0.5;

    Chi_low  = Chi_all(idx_low);
    Chi_high = Chi_all(idx_high);

    figure("Color", "white", "Position", [120, 120, 900, 600]); hold on;

    % Consistent bins
    binWidth = 1;
    if isempty(Chi_all)
        edges = -90:binWidth:90;
    else
        edges = floor(min(Chi_all)):binWidth:ceil(max(Chi_all));
    end
    x_fit = linspace(edges(1), edges(end), 400);

    % Colors
    col_low_data  = [0.25 0.55 0.95];
    col_low_fit   = [0.05 0.25 0.75];
    col_high_data = [0.95 0.55 0.25];
    col_high_fit  = [0.75 0.25 0.05];

    legend_entries = {};
    plot_handles = [];

    % Low group: CoM_y < 0.5
    if ~isempty(Chi_low)
        hL = histogram(Chi_low, 'Normalization', 'pdf', 'BinEdges', edges, ...
            'FaceAlpha', 0.5, 'FaceColor', col_low_data, 'EdgeColor', 'none');
        plot_handles(end+1) = hL;
        legend_entries{end+1} = 'Lower half: $\tilde{y} < 0.5$';

        ChiL = Chi_low(~isnan(Chi_low));
        if isempty(ChiL)
            muL = NaN; sigmaL = NaN; yL = [];
        else
            if numel(ChiL) >= 5
                mu0 = median(ChiL);
                s0  = 1.4826*mad(ChiL, 1);
                if s0 <= eps, s0 = max(std(ChiL,1), eps); end
                k = 2.5;
                inliers = abs(ChiL - mu0) <= k*s0;
                if nnz(inliers) < max(5, round(0.5*numel(ChiL)))
                    inliers = abs(ChiL - mu0) <= 3.5*s0;
                end
                Chi_fitL = ChiL(inliers);
            else
                Chi_fitL = ChiL;
            end
            muL = mean(Chi_fitL);
            sigmaL = std(Chi_fitL, 1);
            if ~isfinite(sigmaL) || sigmaL <= 0, sigmaL = max(std(ChiL,1), 1e-6); end
            yL = (1./(sigmaL*sqrt(2*pi))) .* exp(-0.5*((x_fit - muL)./sigmaL).^2);
        end
        if ~isempty(ChiL)
            pL = plot(x_fit, yL, '-', 'Color', col_low_fit, 'LineWidth', 2);
            plot_handles(end+1) = pL;
            legend_entries{end+1} = sprintf('Fit: $\\mu=%.2f^{\\circ},\\ \\sigma=%.2f^{\\circ}$', muL, sigmaL);
        end
    end

    % High group: CoM_y >= 0.5
    if ~isempty(Chi_high)
        hH = histogram(Chi_high, 'Normalization', 'pdf', 'BinEdges', edges, ...
            'FaceAlpha', 0.5, 'FaceColor', col_high_data, 'EdgeColor', 'none');
        plot_handles(end+1) = hH;
        legend_entries{end+1} = 'Upper half: $\tilde{y} \geq 0.5$';

        ChiH = Chi_high(~isnan(Chi_high));
        if isempty(ChiH)
            muH = NaN; sigmaH = NaN; yH = [];
        else
            if numel(ChiH) >= 5
                mu0 = median(ChiH);
                s0  = 1.4826*mad(ChiH, 1);
                if s0 <= eps, s0 = max(std(ChiH,1), eps); end
                k = 2.5;
                inliers = abs(ChiH - mu0) <= k*s0;
                if nnz(inliers) < max(5, round(0.5*numel(ChiH)))
                    inliers = abs(ChiH - mu0) <= 3.5*s0;
                end
                Chi_fitH = ChiH(inliers);
            else
                Chi_fitH = ChiH;
            end
            muH = mean(Chi_fitH);
            sigmaH = std(Chi_fitH, 1);
            if ~isfinite(sigmaH) || sigmaH <= 0, sigmaH = max(std(ChiH,1), 1e-6); end
            yH = (1./(sigmaH*sqrt(2*pi))) .* exp(-0.5*((x_fit - muH)./sigmaH).^2);
        end
        if ~isempty(ChiH)
            pH = plot(x_fit, yH, '-', 'Color', col_high_fit, 'LineWidth', 2);
            plot_handles(end+1) = pH;
            legend_entries{end+1} = sprintf('Fit: $\\mu=%.2f^{\\circ},\\ \\sigma=%.2f^{\\circ}$', muH, sigmaH);
        end
    end

    xlabel('Orientation $\theta$ (deg)', 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    ylabel('Probability Density', 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    title(sprintf('Actin concentration: %.2f\\%%, Z: %.1f um\n', ...
                   susp_conc, Z_pos), 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    grid on;
    xlim([-25 25]);
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
    legend(plot_handles, legend_entries, 'Location', 'best', 'Interpreter', 'latex', 'FontName', 'Times New Roman');

    % Save the figure
    saveas(gcf, fullfile(exp2proc, sprintf('Chi_PDF_Split_Conc%.2f_Z%.1f.png', susp_conc, Z_pos)));

end


%% Plot S vs Z for different suspension concentrations
% Extract arrays for plotting
Z = [Combined_info.Z]';
S = [Combined_info.S]';
susp_conc = [Combined_info.susp_conc]';

% Plot S vs Z for different Q values
figure("Color", "white", "Position", [100, 100, 800, 600]);
unique_susp_conc = unique(susp_conc);
hold on;
colors = lines(length(unique_susp_conc));
legends = cell(1, length(unique_susp_conc));
for i = 1:length(unique_susp_conc)
    idx = susp_conc == unique_susp_conc(i);
    plot(Z(idx)/500, S(idx), 'o', 'Color', colors(i,:), 'LineWidth', 8, 'MarkerFaceColor', colors(i,:));
    legends{i} = ['Actin concentration = ' num2str(unique_susp_conc(i)) ' %'];
end
hold off;
xlabel('Normalized Z-Position', 'FontName', 'Times New Roman');
ylabel('$S = \frac{1}{2} \left( 3 \langle \cos^2\theta \rangle - 1 \right)$', 'Interpreter', 'latex');
title('Order Parameter S vs Z', 'Interpreter', 'latex');
legend(legends, 'Location', 'southeast');
set(legend, 'FontName', 'Times New Roman');
grid on;
set(gca, 'FontSize', 18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
xlim([0 0.5])

% Save the figure
saveas(gcf, fullfile(exp2proc, 'OrderParameter_vs_Z.png'));


%% Plot average orientation angle for bottom and upper half vs Z and vs concentration
% This section computes axial means (handles 180° symmetry) and circular
% std/SEM for lower (CoM_y < 0.5) and upper (CoM_y >= 0.5) halves, then
% plots mean orientation vs Z for each concentration and mean orientation
% vs concentration for each Z.

% Preallocate arrays
G = numel(Combined_info);
Z_vals = nan(G,1);
conc_vals = nan(G,1);
mean_low = nan(G,1);
mean_high = nan(G,1);
std_low_deg = nan(G,1);
std_high_deg = nan(G,1);
n_low = nan(G,1);
n_high = nan(G,1);

% Helper inline for axial mean/std (degrees)
axial_stats = @(angles_deg) deal( ...
    rad2deg( 0.5 * atan2( mean( sin( 2*deg2rad(angles_deg) ) ), ...
                          mean( cos( 2*deg2rad(angles_deg) ) ) ) ), ... % axial mean (deg)
    rad2deg( 0.5 * sqrt( -2 * log( ...
        max( min( hypot( mean( cos( 2*deg2rad(angles_deg) ) ), ...
                          mean( sin( 2*deg2rad(angles_deg) ) ) ), 1 ), eps ) ) ) ) ... % axial std (deg)
);

for k = 1:G
    Z_vals(k) = Combined_info(k).Z;
    conc_vals(k) = Combined_info(k).susp_conc;
    Chi_all = Combined_info(k).Chi_all(:);
    ynorm = Combined_info(k).CoM_y_all_normalized(:);
    if isempty(Chi_all)
        continue
    end

    idxL = ~isnan(Chi_all) & (ynorm < 0.5);
    idxH = ~isnan(Chi_all) & (ynorm >= 0.5);

    % Lower half
    if any(idxL)
        [mL, sL] = axial_stats(Chi_all(idxL));
        mean_low(k) = mL;
        std_low_deg(k) = sL;
        n_low(k) = nnz(idxL);
    end
    % Upper half
    if any(idxH)
        [mH, sH] = axial_stats(Chi_all(idxH));
        mean_high(k) = mH;
        std_high_deg(k) = sH;
        n_high(k) = nnz(idxH);
    end
end

% Compute SEMs
sem_low = std_low_deg ./ sqrt(max(1,n_low));
sem_high = std_high_deg ./ sqrt(max(1,n_high));

% Unique concentrations and Z positions
unique_conc = unique(conc_vals(~isnan(conc_vals)));
unique_Z = unique(Z_vals(~isnan(Z_vals)));

% 1) Plot mean orientation vs Z for each concentration (lower & upper)
figure('Color','white','Position',[200 200 900 600]); hold on;
colors = lines(numel(unique_conc));
for ic = 1:numel(unique_conc)
    conc = unique_conc(ic);
    idx = conc_vals == conc;
    if ~any(idx), continue; end
    % Sort by Z for nicer lines
    [Zs, sidx] = sort(Z_vals(idx));
    Zs_norm = Zs / 500; % consistent normalization used earlier
    ml = mean_low(idx); ml = ml(sidx);
    mh = mean_high(idx); mh = mh(sidx);
    seml = sem_low(idx); seml = seml(sidx);
    semh = sem_high(idx); semh = semh(sidx);

    % Plot lower half
    errorbar(Zs_norm, ml, seml, 'o-', 'Color', colors(ic,:), 'MarkerFaceColor', colors(ic,:), 'LineWidth', 1.5);
    % Plot upper half with dashed line and different marker
    errorbar(Zs_norm, mh, semh, 's--', 'Color', colors(ic,:)*0.6 + 0.2, 'MarkerFaceColor', colors(ic,:)*0.6 + 0.2, 'LineWidth', 1.5);

    legend_entries{ic*2-1} = sprintf('Conc %.2f%% (lower)', conc);
    legend_entries{ic*2}   = sprintf('Conc %.2f%% (upper)', conc);
end
xlabel('Normalized Z-Position', 'FontName','Times New Roman');
ylabel('Mean Orientation (deg)', 'FontName','Times New Roman');
title('Mean Orientation vs Z (lower vs upper halves)', 'FontName','Times New Roman');
grid on;
set(gca,'FontSize',14,'LineWidth',1.2);
xlim([min(Z_vals(~isnan(Z_vals)))/500, max(Z_vals(~isnan(Z_vals)))/500]);
legend(legend_entries,'Interpreter','none','Location','bestoutside');
hold off;
saveas(gcf, fullfile(exp2proc, 'MeanOrientation_vs_Z_by_Concentration.png'));

% 2) Plot mean orientation vs concentration for each Z (lower & upper)
figure('Color','white','Position',[200 200 900 600]); hold on;
colorsZ = lines(numel(unique_Z));
legend_entries = {};
for iz = 1:numel(unique_Z)
    Zv = unique_Z(iz);
    Zv_norm = Zv / 500;                    % normalized Z
    idx = Z_vals == Zv;
    if ~any(idx), continue; end
    % Sort by concentration for nicer lines
    [Cs, sidx] = sort(conc_vals(idx));
    ml = mean_low(idx); ml = ml(sidx);
    mh = mean_high(idx); mh = mh(sidx);
    seml = sem_low(idx); seml = seml(sidx);
    semh = sem_high(idx); semh = semh(sidx);

    % Lower half
    errorbar(Cs, ml, seml, 'o-', 'Color', colorsZ(iz,:), 'MarkerFaceColor', colorsZ(iz,:), 'LineWidth', 1.5);
    % Upper half
    errorbar(Cs, mh, semh, 's--', 'Color', colorsZ(iz,:)*0.6 + 0.2, 'MarkerFaceColor', colorsZ(iz,:)*0.6 + 0.2, 'LineWidth', 1.5);

    legend_entries{end+1} = sprintf('Z=%.3f (lower)', Zv_norm);
    legend_entries{end+1} = sprintf('Z=%.3f (upper)', Zv_norm);
end
xlabel('Actin concentration (%)', 'FontName','Times New Roman');
ylabel('Mean Orientation (deg)', 'FontName','Times New Roman');
title('Mean Orientation vs Concentration (by Z)', 'FontName','Times New Roman');
grid on;
set(gca,'FontSize',14,'LineWidth',1.2);
legend(legend_entries,'Interpreter','none','Location','bestoutside');
hold off;
saveas(gcf, fullfile(exp2proc, 'MeanOrientation_vs_Concentration_by_Z.png'));

% 3) Optional: table summary
Summary = table(conc_vals, Z_vals, n_low, mean_low, std_low_deg, sem_low, n_high, mean_high, std_high_deg, sem_high, ...
    'VariableNames', {'Concentration','Z','N_low','Mean_low_deg','Std_low_deg','SEM_low_deg','N_high','Mean_high_deg','Std_high_deg','SEM_high_deg'});
save(fullfile(exp2proc, 'MeanOrientation_Summary.mat'), 'Summary');
writetable(Summary, fullfile(exp2proc, 'MeanOrientation_Summary.csv'));