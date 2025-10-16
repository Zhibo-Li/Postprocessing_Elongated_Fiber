%%% actin in PAs: DLD check

clear; close all; clc;

image_height = 204.8; % the height of the image, unit: um
magni = 0.1; % the magnification of the objective (um/pixel)
flowfocusing_bandwidth = 100; % the width of suspension band after flow-focusing, unit: um
X_base = 6000; % unit: um
lambda = [30, 45]; % unit: um

set(0, 'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

fiber_lengths = []; fiber_lengths_normalized = []; fiber_Ys = []; fiber_Xs = [];
for exp_date = 1: 2

    switch exp_date
        case 1
            exp2proc = dir(['E:\Processing & Results\Actin Filaments in Porous Media\' ...
                '20241031-Actin-DLD\Y*']);
            % load extracted data from *.csv files
            longfibers = readmatrix(['E:\Processing & Results\' ...
                'Actin Filaments in Porous Media\20241031-Actin-DLD\' ...
                'Extracted_data_Y-positon_Length_exp20231030_exp20240415_2.csv']);
        case 2
            exp2proc = dir(['E:\Processing & Results\Actin Filaments in Porous Media\' ...
                '20241111-Actin-DLD\Y*']);

    end

    for ii = 1:length(exp2proc)

        Y_base = exp2proc(ii).name;
        Y_base = str2double(cell2mat((extractBetween(Y_base,'Y','um'))));

        Tracings = dir(fullfile(exp2proc(ii).folder, exp2proc(ii).name, 'Tracings_*.csv'));
        Vertices = dir(fullfile(exp2proc(ii).folder, exp2proc(ii).name, 'Vertices_*.csv'));

        for jj = 1:length(Tracings)

            Tracing_info = readmatrix(fullfile(Tracings(jj).folder, Tracings(jj).name),"NumHeaderLines",1);
            Vertice_info = readcell(fullfile(Vertices(jj).folder, Vertices(jj).name),"NumHeaderLines",1);

            fiber_lengths = [fiber_lengths; Tracing_info(:, 6)*magni];
            fiber_index = Vertice_info(:,2);
            fiber_index = cell2mat(cellfun(@(x) str2double(x(2:end)), fiber_index, 'UniformOutput', false));

            for kk = 1:max(fiber_index)

                fiber_xy = cell2mat(Vertice_info(fiber_index == kk,5:6));
                fiber_CoM_xy = mean(fiber_xy,1)*magni;

                fiber_Ys = [fiber_Ys; image_height-fiber_CoM_xy(2)+Y_base];
                fiber_Xs = [fiber_Xs; fiber_CoM_xy(1)+X_base];

            end
        end
    end

    if exp_date == 1
        fiber_Xs = [fiber_Xs; ones(size(longfibers(:,2)))*9000];
        fiber_Ys = [fiber_Ys; longfibers(:,2)];
        fiber_lengths = [fiber_lengths; longfibers(:,1)];

        data_num = size(longfibers(:,1), 1);
    end

    fiber_lengths_normalized = fiber_lengths/lambda(exp_date);

    % calculate the drift angle of the fibers
    fiber_Ys_max = fiber_Ys;
    fiber_Ys_min = fiber_Ys - flowfocusing_bandwidth;
    if exp_date == 1
        fiber_Ys_max(end-data_num+1:end) = fiber_Ys(end-data_num+1:end) + 200;
        fiber_Ys_min(end-data_num+1:end) = fiber_Ys(end-data_num+1:end) - 200;
    end
    fiber_Ys_min(fiber_Ys_min < 0) = 0; % set value to zeros if it's negative
    beta_max{exp_date} = atand(fiber_Ys_max./fiber_Xs);
    beta_min{exp_date} = atand(fiber_Ys_min./fiber_Xs);
    fiber_lengths_together{exp_date} = fiber_lengths;
    fiber_lengths_together_normalized{exp_date} = fiber_lengths_normalized;

    % bin average for errorbar plot (based on fiber length)
    n = 1;
    for ii = 0:lambda(exp_date)/5:100

        if ii < 80
            YY_max = fiber_Ys_max(fiber_lengths > ii & fiber_lengths < ii+10);
            YY_min = fiber_Ys_min(fiber_lengths > ii & fiber_lengths < ii+10);
            XX = fiber_Xs(fiber_lengths > ii & fiber_lengths < ii+10);
            LL = fiber_lengths(fiber_lengths > ii & fiber_lengths < ii+10);
        else
            YY_max = fiber_Ys_max(fiber_lengths > ii);
            YY_min = fiber_Ys_min(fiber_lengths > ii);
            XX = fiber_Xs(fiber_lengths > ii);
            LL = fiber_lengths(fiber_lengths > ii);
        end
        Y_max_mean(exp_date,n) = mean(YY_max); Y_max_std(exp_date,n) = std(YY_max);
        Y_min_mean(exp_date,n) = mean(YY_min); Y_min_std(exp_date,n) = std(YY_min);
        X_mean(exp_date,n) = mean(XX); X_std(exp_date,n) = std(XX);
        L_mean(exp_date,n) = mean(LL); L_std(exp_date,n) = std(LL);
        n=n+1;

    end

    % bin average for errorbar plot (based on normalized fiber length)
    n = 1; bin_size = 0.12; bin_number = 20;
    for ii = 0:bin_number

        if ii < bin_number
            YY_max = fiber_Ys_max(fiber_lengths_normalized > ii*bin_size & fiber_lengths_normalized < (ii+1)*bin_size);
            YY_min = fiber_Ys_min(fiber_lengths_normalized > ii*bin_size & fiber_lengths_normalized < (ii+1)*bin_size);
            XX = fiber_Xs(fiber_lengths_normalized > ii*bin_size & fiber_lengths_normalized < (ii+1)*bin_size);
            LL = fiber_lengths_normalized(fiber_lengths_normalized > ii*bin_size & fiber_lengths_normalized < (ii+1)*bin_size);
        else
            YY_max = fiber_Ys_max(fiber_lengths_normalized > ii*bin_size);
            YY_min = fiber_Ys_min(fiber_lengths_normalized > ii*bin_size);
            XX = fiber_Xs(fiber_lengths_normalized > ii*bin_size);
            LL = fiber_lengths_normalized(fiber_lengths_normalized > ii*bin_size);
        end
        Y_max_mean_normalized(exp_date,n) = mean(YY_max); Y_max_std_normalized(exp_date,n) = std(YY_max);
        Y_min_mean_normalized(exp_date,n) = mean(YY_min); Y_min_std_normalized(exp_date,n) = std(YY_min);
        X_mean_normalized(exp_date,n) = mean(XX); X_std_normalized(exp_date,n) = std(XX);
        L_mean_normalized(exp_date,n) = mean(LL); L_std_normalized(exp_date,n) = std(LL);
        n=n+1;

    end

    % initialize for next loop
    fiber_lengths = []; fiber_lengths_normalized = [];  fiber_Ys = []; fiber_Xs = [];
end

figure('Color','w','Position',[100 100 800 600]);
% first plot: exp20241031
plot(fiber_lengths_together_normalized{1}, beta_max{1}, 'color', [203, 213, 232]/255, ...
    'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 8, 'LineWidth', 1.2);
hold on;
% second plot: exp20241111
plot(fiber_lengths_together_normalized{2}, beta_max{2}, 'color', [253, 205, 172]/255, ...
    'LineStyle', 'none', 'Marker', '^', 'MarkerSize', 8, 'LineWidth', 1.2);
hold on;

% % % %% errorbar plot: bin average based on fiber length
% % % % errorbar plot: exp20241031
% % % errorbar(L_mean(1, :)/lambda(1), atand(Y_max_mean(1, :)./X_mean(1, :)), atand(Y_max_std(1, :)./X_mean(1, :)), ...
% % %     atand(Y_max_std(1, :)./X_mean(1, :)), L_std(1, :)/lambda(1), L_std(1, :)/lambda(1), 'color', [27, 158, 119]/255, ...
% % %     'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2);
% % % hold on;
% % % % errorbar plot: exp20241111
% % % errorbar(L_mean(2, :)/lambda(2), atand(Y_max_mean(2, :)./X_mean(2, :)), atand(Y_max_std(2, :)./X_mean(2, :)), ...
% % %     atand(Y_max_std(2, :)./X_mean(2, :)), L_std(2, :)/lambda(2), L_std(2, :)/lambda(2), 'color', [217, 95, 2]/255, ...
% % %     'LineStyle', 'none', 'Marker', '^', 'MarkerSize', 10, 'LineWidth', 2);
% % % hold on;

%% errorbar plot: bin average based on fiber length normalized
% errorbar plot: exp20241031
errorbar(L_mean_normalized(1, :), atand(Y_max_mean_normalized(1, :)./X_mean_normalized(1, :)), atand(Y_max_std_normalized(1, :)./X_mean_normalized(1, :)), ...
    atand(Y_max_std_normalized(1, :)./X_mean_normalized(1, :)), L_std_normalized(1, :), L_std_normalized(1, :), 'color', [117, 112, 179]/255, ...
    'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
% errorbar plot: exp20241111
errorbar(L_mean_normalized(2, :), atand(Y_max_mean_normalized(2, :)./X_mean_normalized(2, :)), atand(Y_max_std_normalized(2, :)./X_mean_normalized(2, :)), ...
    atand(Y_max_std_normalized(2, :)./X_mean_normalized(2, :)), L_std_normalized(2, :), L_std_normalized(2, :), 'color', [217, 95, 2]/255, ...
    'LineStyle', 'none', 'Marker', '^', 'MarkerSize', 10, 'LineWidth', 2);
hold on;

ylabel('$\beta\ (\mathrm{^\circ})$','FontSize', 24,'Interpreter', 'latex');
xlabel('$L/\lambda$','FontSize', 24,'Interpreter', 'latex');
xlim([0 3.3]); ylim([0 10])
set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24)

legend({'', '', '$\lambda=30\mathrm{\mu m}$', ...
    '$\lambda=45\mathrm{\mu m}$'},'Interpreter', 'latex', 'FontSize', 24)

text(2.6, 0.7, '$\lambda/R=3$', 'FontSize', 24, 'FontWeight', 'bold', ...
    'Interpreter','latex','BackgroundColor','w');

f=gcf;
saveas(f, ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Actin-DLD\Beta_Length_errorbarPlot_exp20241031_exp20241111_BinOnNormalized.fig'])

f=gcf;
exportgraphics(f,['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Actin-DLD\Beta_Length_errorbarPlot_exp20241031_exp20241111_BinOnNormalized.png'],'Resolution',100)

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Actin-DLD\Beta_Length_errorbarPlot_exp20241031_exp20241111_BinOnNormalized.eps']);


