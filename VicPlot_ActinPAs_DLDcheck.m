%%% actin in PAs: DLD check

clear; close all; clc;

image_height = 204.8; % the height of the image, unit: um 
magni = 0.1; % the magnification of the objective (um/pixel)

set(0, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

fiber_lengths = []; fiber_Ys = [];
for exp_date = 1:2 % 3

    switch exp_date
        case 1
            exp2proc = dir(['F:\Processing & Results\Actin Filaments in Porous Media\' ...
                '20231030-Actin-DLD\Y*']);
        case 2
            exp2proc = dir(['F:\Processing & Results\Actin Filaments in Porous Media\' ...
                '20240415-Actin-DLD\Y*']);
        case 3
            exp2proc = dir(['F:\Processing & Results\Actin Filaments in Porous Media\' ...
                '20240522-Actin-DLD\Y*']);
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

            end
        end
    end
end

%% plot 1
figure('Color','w','Position',[100 100 800 600]); 
plot(fiber_lengths,fiber_Ys, 'r.', 'MarkerSize', 28);
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


%% errorbar plot
n = 1;
for ii = 5:10:100

    if ii < 80
        YY = fiber_Ys(fiber_lengths > ii & fiber_lengths < ii+10);
        LL = fiber_lengths(fiber_lengths > ii & fiber_lengths < ii+10);
    else
        YY = fiber_Ys(fiber_lengths > ii);
        LL = fiber_lengths(fiber_lengths > ii);
    end
    L_mean(n) = mean(LL); L_std(n) = std(LL);
    Y_mean(n) = mean(YY); Y_std(n) = std(YY);
    n=n+1;

end

figure('Color','w','Position',[100 100 800 600]); 
plot(fiber_lengths,fiber_Ys, 'm.', 'MarkerSize', 12);
hold on;
errorbar(L_mean, Y_mean, Y_std, Y_std, L_std, L_std, 'ok', 'LineStyle', ...
    'none', 'LineWidth',2);
hold on;
ylabel('$y\ (\mathrm{\mu m})$','FontSize', 24,'Interpreter', 'latex');
xlabel('$L\ (\mathrm{\mu m})$','FontSize', 24,'Interpreter', 'latex');
xlim([0 100]); ylim([0 1500])
set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24)

f=gcf;
exportgraphics(f,['F:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Actin-DLD\Y-positon_Length_errorbarPlot_exp20231030_exp20240415.png'],'Resolution',100)

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Actin-DLD\Y-positon_Length_errorbarPlot_exp20231030_exp20240415.eps']);
