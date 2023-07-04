%%% Plotting for actin filament vs triangular pillar

clear; close all; clc;

% basic information.
Obj_Mag = 0.1;

% laod data and put all in one.
xlsfile_1 = readcell(['F:\Processing & Results\FSI - Actin &  Individual Obstacle\' ...
    'Exp data triangular pillar lateral pointing\20230330-Actin-Individual_' ...
    'triangularPillar_uppoint.xlsx'],'Sheet','Sheet1','NumHeaderLines',1);
the_names_1 = xlsfile_1(:, 1);
the_data_1 = xlsfile_1(:, 2:6);
mask = cellfun(@ismissing, the_data_1); the_data_1(mask) = {nan};
the_data_1 = cell2mat(the_data_1); 

xlsfile_2 = readcell(['F:\Processing & Results\FSI - Actin &  Individual Obstacle\' ...
    'Exp data triangular pillar lateral pointing\20230405-Actin-Individual_' ...
    'triangularPillar_uppoint.xlsx'],'Sheet','Sheet1','NumHeaderLines',1);
the_names_2 = xlsfile_2(:, 1);
the_data_2 = xlsfile_2(:, 2:6);
mask = cellfun(@ismissing, the_data_2); the_data_2(mask) = {nan};
the_data_2 = cell2mat(the_data_2);

the_names = [the_names_1; the_names_2];
the_data = [the_data_1; the_data_2];

% %     No.1 column: the_L (um)
% %     No.2 column: normalized the_y0
% %     No.3 column: normalized the_delta0
% %     No.4 column: the_theta0
% %     No.5 column: elastoviscous number

% % % % %% statistics
% % % % % Contour length L
% % % % figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% % % % histogram(the_data(:,1));
% % % % set(gca,'FontSize',16);
% % % % xlabel('$Contour\ length\ L$','FontSize', 22,'Interpreter', 'latex');
% % % % ylabel('$Number\ of\ cases$','FontSize', 22,'Interpreter', 'latex');
% % % % % xlim([0 2.2]);
% % % % % f=gcf;
% % % % % exportgraphics(f,'Statistics_contourL.png','Resolution',100)
% % % % 
% % % % % Initial angle chi_0
% % % % figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% % % % edges = -90:20:90;
% % % % histogram(the_data(:,4),edges);
% % % % set(gca,'FontSize',16);
% % % % xlabel('$Initial\ angle\ \theta_0\ (^{\circ})$','FontSize', 22,'Interpreter', 'latex');
% % % % ylabel('$Number\ of\ cases$','FontSize', 22,'Interpreter', 'latex');
% % % % xlim([-90 90]);
% % % % % f=gcf;
% % % % % exportgraphics(f,'Statistics_Chi0.png','Resolution',100)
% % % % 
% % % % % Normalized initial position (y_0/h_obs)
% % % % figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% % % % edges = -0.2:0.2:1.2;
% % % % histogram(the_data(:,2),edges);
% % % % set(gca,'FontSize',16);
% % % % xlabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% % % % ylabel('$Number\ of\ cases$','FontSize', 22,'Interpreter', 'latex');
% % % % xlim([-0.2 1.2]);
% % % % % f=gcf;
% % % % % exportgraphics(f,'Statistics_y0.png','Resolution',100)

% create dialog box to gather user input for plotting;
names_plot = the_names'; data_plot = the_data';
prompt = {'The lower bound of the initial angle:', 'The upper bound of the initial angle:', ...
    'The lower bound of the contour length:','The upper bound of the contour length:'...
    'The lower bound of the initial position:','The upper bound of the initial position:'};
definput = {'nan', 'nan', 'nan', 'nan', '0', '1'};
% definput = {'nan', 'nan', 'nan', 'nan', 'nan', 'nan'};
answer = inputdlg(prompt, 'Input (please input NaN if there is no bound)', [1 35] , definput);

%%%%%%%%%%%%%%%%%%%%%%%%% assign the values %%%%%%%%%%%%%%%%%%%%%%%%%%
range_chi0_low = str2double(answer{1,1}); if isnan(range_chi0_low); range_chi0_low = -91; end
range_chi0_up = str2double(answer{2,1}); if isnan(range_chi0_up); range_chi0_up = 91; end
range_L_low = str2double(answer{3,1}); if isnan(range_L_low); range_L_low = 0; end
range_L_up = str2double(answer{4,1}); if isnan(range_L_up); range_L_up = 100; end
range_y0_low = str2double(answer{5,1}); if isnan(range_y0_low); range_y0_low = -10; end
range_y0_up = str2double(answer{6,1}); if isnan(range_y0_up); range_y0_up = 10; end
% chi_0
names_plot(data_plot(4, :) < range_chi0_low) = [];
data_plot(:, data_plot(4, :) < range_chi0_low) = []; 
names_plot(data_plot(4, :) > range_chi0_up) = [];
data_plot(:, data_plot(4, :) > range_chi0_up) = [];
% L
names_plot(data_plot(1, :) < range_L_low) = [];
data_plot(:, data_plot(1, :) < range_L_low) = []; 
names_plot(data_plot(1, :) > range_L_up) = [];
data_plot(:, data_plot(1, :) > range_L_up) = [];
% y_0
names_plot(data_plot(2, :) < range_y0_low) = [];
data_plot(:, data_plot(2, :) < range_y0_low) = [];  
names_plot(data_plot(2, :) > range_y0_up) = [];
data_plot(:, data_plot(2, :) > range_y0_up) = [];

data_plot(1, :) = data_plot(1, :)/25; % normalize contour length by l_obs (25 um)

%% plotting
%%%%%%%%%%%%%%%%%%%%%%%%% scatter plots %%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the deviation vs L (mubar) & y_0 
figure('color', 'w'); 
set(gcf, 'Position',[100 100 1000 600], 'Color','white', 'DefaultTextInterpreter', 'latex')

scatter(data_plot(2, :)', data_plot(1, :)', 200, data_plot(3, :)', 'Filled','MarkerEdgeColor','k'); hold on
set(gca, 'YScale', 'log', 'YLim', [0 2], 'YTick', [0.25 0.5 1 2], 'YTickLabel', {'$0.25$','$0.5$', '$1$', '$2$'})
xlabel('$y_0/h_{\rm obs}$','FontSize', 24,'Interpreter', 'latex');
ylabel('$L/l_{\rm obs}$','FontSize', 24,'Interpreter', 'latex');

caxis([-0.1 0.1]); cmocean('balance')
hcb=colorbar; 
t = title(hcb,'$\delta/h_{\rm obs}$','FontSize', 24,'Interpreter', 'latex');
% hcb.Label.String = '$\delta/h_{\rm obs}$';
% hcb.Label.Interpreter = 'LaTeX';
hcb.TickLabelInterpreter = 'LaTeX';
% hcb.FontSize = 24;
hcb.Position = [0.9, 0.15, 0.015, 0.66];
t.Position = [11.6,327,0];

yyaxis right
plot(data_plot(2, :)', data_plot(5, :)', 'LineStyle','none');
set(gca, 'YScale', 'log', 'YLim', [0  VicFc_Get_elastoviscousNum(50*1e-6, 2e-4, 2.5e-5)]) % YLim should be aligned with the length
ylabel('$\bar{\mu}$','FontSize', 24,'Interpreter', 'latex');

set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24, ...
    'NextPlot','replacechildren', 'TickLabelInterpreter','latex', 'Position', [0.13,0.15,0.61,0.8])

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
    'FSI - Actin &  Individual Obstacle\Figure\deviation_y0_L_mubar.eps']);

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['F:\Processing & Results\FSI - Actin &  Individual Obstacle' ...
    '\Figure\deviation_y0_L_mubar.pdf']);
