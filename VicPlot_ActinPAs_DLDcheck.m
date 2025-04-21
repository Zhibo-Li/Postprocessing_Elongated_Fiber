%%% actin in PAs: DLD check

clear; close all; clc;

image_height = 204.8; % the height of the image, unit: um 
magni = 0.1; % the magnification of the objective (um/pixel)
flowfocusing_bandwidth = 100; % the width of suspension band after flow-focusing, unit: um

set(0, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

fiber_lengths = []; fiber_Ys = []; fiber_Xs = [];
for exp_date =  5 % 1:2 % 3

    if exp_date == 1 || exp_date == 2 || exp_date == 4
        lambda = 30; % center-to-center distance of pillars in pillar array
    else 
        lambda = 45; % center-to-center distance of pillars in pillar array
    end

    switch exp_date
        case 1 % data lost
            exp2proc = dir(['F:\Processing & Results\Actin Filaments in Porous Media\' ...
                '20231030-Actin-DLD\Y*']);
        case 2 % data lost
            exp2proc = dir(['F:\Processing & Results\Actin Filaments in Porous Media\' ...
                '20240415-Actin-DLD\Y*']);
        case 3 % data lost
            exp2proc = dir(['F:\Processing & Results\Actin Filaments in Porous Media\' ...
                '20240522-Actin-DLD\Y*']);
        case 4
            exp2proc = dir(['E:\Processing & Results\Actin Filaments in Porous Media\' ...
                '20241031-Actin-DLD\Y*']);
            X_base = 6000; % unit: um
            % load extracted data from *.csv files
            longfibers = readmatrix(['E:\Processing & Results\' ...
                'Actin Filaments in Porous Media\20241031-Actin-DLD\' ...
                'Extracted_data_Y-positon_Length_exp20231030_exp20240415.csv']);
        case 5
            exp2proc = dir(['E:\Processing & Results\Actin Filaments in Porous Media\' ...
                '20241111-Actin-DLD\Y*']);
            X_base = 6000; % unit: um
            
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

    if exp_date == 4
        fiber_Xs = [fiber_Xs; ones(size(longfibers(:,2)))*9000];
        fiber_Ys = [fiber_Ys; longfibers(:,2)];
        fiber_lengths = [fiber_lengths; longfibers(:,1)];

        data_num = size(longfibers(:,1), 1);
    end

end

%% plot 1
figure('Color','w','Position',[100 100 800 600]); 
plot(fiber_lengths/45,fiber_Ys, 'r.', 'MarkerSize', 28);
ylabel('$y\ (\mathrm{\mu m})$','FontSize', 24,'Interpreter', 'latex');
xlabel('$L\ (\mathrm{\mu m})$','FontSize', 24,'Interpreter', 'latex');
set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24)

f=gcf;
exportgraphics(f,['F:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Actin-DLD\Y-positon_Length.png'],'Resolution',100)


%% violin plot
n = 1;
for ii = 5:10:55

    if ii < 50
        YY{n} = fiber_Ys(fiber_lengths > ii & fiber_lengths < ii+10);
        LL = fiber_lengths(fiber_lengths > ii & fiber_lengths < ii+10);
    else
        YY{n} = fiber_Ys(fiber_lengths > ii);
        LL = fiber_lengths(fiber_lengths > ii);
    end
    L_mean(n) = mean(LL); L_std(n) = std(LL);
    n=n+1;

end

figure('Color','w','Position',[100 100 800 600]); 
errorbar(10:10:60, 1250*ones(1, 6),L_std,'ok', 'horizontal', 'LineStyle', ...
    'none', 'LineWidth',2);
hold on;
plot(fiber_lengths,fiber_Ys, 'm.', 'MarkerSize', 12);
hold on;
violin(YY,'x',[10 20 30 40 50 60]); 
ylabel('$y\ (\mathrm{\mu m})$','FontSize', 24,'Interpreter', 'latex');
xlabel('$L\ (\mathrm{\mu m})$','FontSize', 24,'Interpreter', 'latex');
xlim([0 70]); ylim([0 1500])
set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24)

f=gcf;
exportgraphics(f,['F:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Actin-DLD\Y-positon_Length_violinPlot.png'],'Resolution',100)

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Actin-DLD\Y-positon_Length_violinPlot.eps']);


%% errorbar plot (bin average)
figure('Color','w','Position',[100 100 800 600]);
fiber_Ys_max = fiber_Ys;
fiber_Ys_min = fiber_Ys - flowfocusing_bandwidth;
fiber_Ys_min(end-data_num+1:end) = fiber_Ys(end-data_num+1:end) - 200; % different band width: 200 um
fiber_Ys_min(fiber_Ys_min < 0) = 0; % set value to zeros if it's negative
beta_max = atand(fiber_Ys_max./fiber_Xs);
beta_min = atand(fiber_Ys_min./fiber_Xs);

plot(fiber_lengths/lambda, beta_max, 'color', [179, 226, 205]/255, ...
    'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 6, 'LineWidth', 1.5);
hold on;

n = 1;
for ii = 5:5:100

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
    Y_max_mean(n) = mean(YY_max); Y_max_std(n) = std(YY_max);
    Y_min_mean(n) = mean(YY_min); Y_min_std(n) = std(YY_min);
    X_mean(n) = mean(XX); X_std(n) = std(XX);
    L_mean(n) = mean(LL); L_std(n) = std(LL);
    n=n+1;

end

errorbar(L_mean/lambda, atand(Y_max_mean./X_mean), atand(Y_max_std./X_mean), ...
    atand(Y_max_std./X_mean), L_std/lambda, L_std/lambda, 'color', [27, 158, 119]/255, ...
        'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
% % % errorbar(L_mean/lambda, atand(Y_min_mean./X_mean), atand(Y_min_std./X_mean), ...
% % %     atand(Y_min_std./X_mean), L_std/lambda, L_std/lambda, 'ok', 'LineStyle', ...
% % %     'none', 'LineWidth',2);
% % % hold on;


%%% errorbar plot: overlap on the previous plot 
fiber_Ys_max = fiber_Ys;
fiber_Ys_min = fiber_Ys - flowfocusing_bandwidth;
fiber_Ys_min(fiber_Ys_min < 0) = 0; % set value to zeros if it's negative
beta_max = atand(fiber_Ys_max./fiber_Xs);
beta_min = atand(fiber_Ys_min./fiber_Xs);

plot(fiber_lengths/lambda, beta_max, 'color', [253, 205, 172]/255, ...
    'LineStyle', 'none', 'Marker', '^', 'MarkerSize', 6, 'LineWidth', 1.5);
hold on;

n = 1;
for ii = 5:5:100

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
    Y_max_mean(n) = mean(YY_max); Y_max_std(n) = std(YY_max);
    Y_min_mean(n) = mean(YY_min); Y_min_std(n) = std(YY_min);
    X_mean(n) = mean(XX); X_std(n) = std(XX);
    L_mean(n) = mean(LL); L_std(n) = std(LL);
    n=n+1;

end

errorbar(L_mean/lambda, atand(Y_max_mean./X_mean), atand(Y_max_std./X_mean), ...
    atand(Y_max_std./X_mean), L_std/lambda, L_std/lambda, 'color', [217, 95, 2]/255, ...
        'LineStyle', 'none', 'Marker', '^', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
% % % errorbar(L_mean/lambda, atand(Y_min_mean./X_mean), atand(Y_min_std./X_mean), ...
% % %     atand(Y_min_std./X_mean), L_std/lambda, L_std/lambda, 'ok', 'LineStyle', ...
% % %     'none', 'LineWidth',2);
% % % hold on;
ylabel('$\beta\ (\mathrm{^\circ})$','FontSize', 24,'Interpreter', 'latex');
xlabel('$L/\lambda$','FontSize', 24,'Interpreter', 'latex');
xlim([0 3.3]); ylim([0 10])
set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24)

legend({'$\lambda=30\mathrm{\mu m}$', ...
    '$\lambda=30\mathrm{\mu m}$', ...
    '$\lambda=45\mathrm{\mu m}$', ...
    '$\lambda=45\mathrm{\mu m}$'},'Interpreter', 'latex', 'FontSize', 18)

text(2.5, 0.7, '$\lambda/a=1.5$', 'FontSize', 24, 'FontWeight', 'bold', ...
    'Interpreter','latex','BackgroundColor','w');

f=gcf;
exportgraphics(f,['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Actin-DLD\Beta_Length_errorbarPlot_exp20241031_exp20241111.png'],'Resolution',100)

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Actin-DLD\Beta_Length_errorbarPlot_exp20241031_exp20241111.eps']);



%% bar plot and errorbar plot
figure('Color','w','Position',[100 100 800 600]);
fiber_Ys_max = fiber_Ys;
fiber_Ys_min = fiber_Ys - flowfocusing_bandwidth;
fiber_Ys_min(end-data_num+1:end) = fiber_Ys(end-data_num+1:end) - 200; % different band width: 200 um
fiber_Ys_min(fiber_Ys_min < 0) = 0; % set value to zeros if it's negative
beta_max = atand(fiber_Ys_max./fiber_Xs);
beta_min = atand(fiber_Ys_min./fiber_Xs);

plot([fiber_lengths'/lambda; fiber_lengths'/lambda], [beta_min'; beta_max'], 'color', [132 186 91]/255, 'LineWidth', 2);
hold on;

n = 1;
for ii = 5:10:100

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
    Y_max_mean(n) = mean(YY_max); Y_max_std(n) = std(YY_max);
    Y_min_mean(n) = mean(YY_min); Y_min_std(n) = std(YY_min);
    X_mean(n) = mean(XX); X_std(n) = std(XX);
    L_mean(n) = mean(LL); L_std(n) = std(LL);
    n=n+1;

end

errorbar(L_mean/lambda, atand(Y_max_mean./X_mean), atand(Y_max_std./X_mean), ...
    atand(Y_max_std./X_mean), L_std/lambda, L_std/lambda, 'ok', 'LineStyle', ...
    'none', 'LineWidth',2);
hold on;
errorbar(L_mean/lambda, atand(Y_min_mean./X_mean), atand(Y_min_std./X_mean), ...
    atand(Y_min_std./X_mean), L_std/lambda, L_std/lambda, 'ok', 'LineStyle', ...
    'none', 'LineWidth',2);
hold on;
% ylabel('$y\ (\mathrm{\mu m})$','FontSize', 24,'Interpreter', 'latex');
% xlabel('$L\ (\mathrm{\mu m})$','FontSize', 24,'Interpreter', 'latex');
ylabel('$\beta\ (\mathrm{^\circ})$','FontSize', 24,'Interpreter', 'latex');
xlabel('$L/\lambda$','FontSize', 24,'Interpreter', 'latex');
xlim([0 3.3]); ylim([0 10])
set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24)

f=gcf;
exportgraphics(f,['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Actin-DLD\Beta_Length_errorbarPlot_exp20241031_plus_VeryLongOnes.png'],'Resolution',100)

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['E:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Actin-DLD\Beta_Length_errorbarPlot_exp20241031_VeryLongOnes.eps']);