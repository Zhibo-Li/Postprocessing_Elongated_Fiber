%%% Plotting for actin filament vs triangular pillar

clear; close all; clc;

% laod data and put all in one.
xlsfile_1 = readcell(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'results_actin_filaments_2023_06_30_Zhibo.xlsx'],'Sheet','Sheet1','NumHeaderLines',1);
the_data = cell2mat(xlsfile_1(:, 1:7));

% %     No.1 column: theta0
% %     No.2 column: normalized y0
% %     No.3 column: normalized L 
% %     No.4 column: elastoviscous number  
% %     No.5 column: delta
% %     No.7 column: elastoviscous number (consider slenderness c)

% create dialog box to gather user input for plotting;
data_plot = the_data';
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
data_plot(:, data_plot(1, :) < range_chi0_low) = []; 
data_plot(:, data_plot(1, :) > range_chi0_up) = [];
% L
data_plot(:, data_plot(3, :) < range_L_low) = []; 
data_plot(:, data_plot(3, :) > range_L_up) = [];
% y_0
data_plot(:, data_plot(2, :) < range_y0_low) = [];  
data_plot(:, data_plot(2, :) > range_y0_up) = [];



%% plotting
%%%%%%%%%%%%%%%%%%%%%%%%% scatter plots %%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the deviation vs L (mubar) & y_0 (Notice the theta_0)
figure;
set(gcf, 'Position',[100 100 1000 600], 'Color','white', 'DefaultTextInterpreter', 'latex')

data_plot_close_to_exp_L0o5_mubar1e4 = data_plot;
data_plot_close_to_exp_L0o5_mubar1e4(:, (data_plot_close_to_exp_L0o5_mubar1e4(4, :) ~= 1e5) ...
    | (data_plot_close_to_exp_L0o5_mubar1e4(3, :) ~= 0.5)) = [];
% choose the numerical data that close to actin property: L = 12um (L/l_obs = 0.5); mu_bar = 3.2e4
% % % data_plot_close_to_exp_L0o75_mubar1e5 = data_plot;
% % % data_plot_close_to_exp_L0o75_mubar1e5(:, (data_plot_close_to_exp_L0o75_mubar1e5(4, :) ~= 1e5) ...
% % %     | (data_plot_close_to_exp_L0o75_mubar1e5(3, :) ~= 0.75)) = [];
% % % % choose the numerical data that close to actin property: : L = 18um (L/l_obs = 0.75); mu_bar = 1.5e5
data_plot_NOT_close_to_exp = data_plot;
data_plot_NOT_close_to_exp(:, (data_plot_NOT_close_to_exp(4, :) == 1e5) ...
    & (data_plot_NOT_close_to_exp(3, :) == 0.5)) = [];
% % % data_plot_NOT_close_to_exp(:, (data_plot_NOT_close_to_exp(4, :) == 1e5) ...
% % %     & (data_plot_NOT_close_to_exp(3, :) == 0.75)) = [];
% the numerical data that NOT close to actin property

scatter3(data_plot_close_to_exp_L0o5_mubar1e4(7, :)', data_plot_close_to_exp_L0o5_mubar1e4(3, :)', ...
    data_plot_close_to_exp_L0o5_mubar1e4(2, :)', 300, data_plot_close_to_exp_L0o5_mubar1e4(5, :)', ...
    'o','Filled','MarkerEdgeColor','k'); hold on
% % % scatter3(data_plot_close_to_exp_L0o75_mubar1e5(7, :)', data_plot_close_to_exp_L0o75_mubar1e5(3, :)', ...
% % %     data_plot_close_to_exp_L0o75_mubar1e5(2, :)', 300, data_plot_close_to_exp_L0o75_mubar1e5(5, :)', ...
% % %     'o','Filled','MarkerEdgeColor','k'); hold on
scatter3(data_plot_NOT_close_to_exp(7, :)', data_plot_NOT_close_to_exp(3, :)', ...
    data_plot_NOT_close_to_exp(2, :)', 300, data_plot_NOT_close_to_exp(5, :)', ...
    'square','Filled','MarkerEdgeColor','k'); hold on
set(gca, 'XScale', 'log', 'XLim', [5 1e4])
xlabel('$\bar{\mu}$','FontSize', 24,'Interpreter', 'latex');
ylabel('$L/l_{\rm obs}$','FontSize', 24,'Interpreter', 'latex');
zlabel('$y_0/h_{\rm obs}$','FontSize', 24,'Interpreter', 'latex');

caxis([-0.1 0.1]); cmocean('balance')
hcb=colorbar; 
hcb.Label.String = '$\delta/h_{\rm obs}$';
hcb.Label.Interpreter = 'LaTeX';
hcb.TickLabelInterpreter = 'LaTeX';
hcb.FontSize = 24;
% hcb.Location = 'eastoutside';
hcb.Position = [0.86, 0.2, 0.015, 0.66];

set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 24, ...
    'TickLabelInterpreter','latex', 'Position', [0.13,0.15,0.63,0.8])

view(-11.71, 4.55)

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
    'FSI - Actin &  Individual Obstacle\Figure\deviation_y0_L_mubar_simu_theta_0.eps']);

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['F:\Processing & Results\FSI - Actin &  Individual Obstacle' ...
    '\Figure\deviation_y0_L_mubar_simu_theta_0.pdf']);
