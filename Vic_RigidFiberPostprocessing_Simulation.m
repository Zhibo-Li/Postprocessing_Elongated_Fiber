% Simulation data plotting.
% Information is stored in the Excel file.



%% # About interaction index 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;
xlsfile = readcell(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\results_2023_02_23_Zhibo.xlsx'],'Sheet','Sheet1' ...
    ,'NumHeaderLines',1);
mask = cellfun(@ismissing, xlsfile); xlsfile(mask) = {nan};
together_plot = [cell2mat(xlsfile(:, 1:5)), cell2mat(xlsfile(:, 12))]; 
together_plot = together_plot';
% thedeg = cell2mat(xlsfile(:, 1)); 
% theL = round(cell2mat(xlsfile(:, 2)), 1); 
% they0 = round(cell2mat(xlsfile(:, 3)), 3); 
% delta = cell2mat(xlsfile(:, 4));
% interaction2 = cell2mat(xlsfile(:, 11)); 
prompt = {'The lower bound of the initial angle:', 'The upper bound of the initial angle:', ...
    'The lower bound of the contour length:','The upper bound of the contour length:'...
    'The lower bound of the initial position:','The upper bound of the initial position:'};
definput = {'-10', '10', 'nan', 'nan', 'nan', 'nan'};
answer = inputdlg(prompt, 'Input (please input NaN if there is no bound)', [1 35] , definput);

%%% assign the values
range_chi0_low = str2double(answer{1,1}); if isnan(range_chi0_low); range_chi0_low = -91; end
range_chi0_up = str2double(answer{2,1}); if isnan(range_chi0_up); range_chi0_up = 91; end
range_L_low = str2double(answer{3,1}); if isnan(range_L_low); range_L_low = 0; end
range_L_up = str2double(answer{4,1}); if isnan(range_L_up); range_L_up = 10; end
range_y0_low = str2double(answer{5,1}); if isnan(range_y0_low); range_y0_low = -10; end
range_y0_up = str2double(answer{6,1}); if isnan(range_y0_up); range_y0_up = 10; end
% chi_0
together_plot(:, together_plot(1, :) < range_chi0_low) = []; 
together_plot(:, together_plot(1, :) > range_chi0_up) = [];
% L
together_plot(:, together_plot(2, :) < range_L_low) = [];
together_plot(:, together_plot(2, :) > range_L_up) = [];
% y_0
together_plot(:, together_plot(3, :) < range_y0_low) = [];  
together_plot(:, together_plot(3, :) > range_y0_up) = [];


figure('color', 'w'); set(gcf, 'Position', [100 100 1000 600]);
% cmap = cmocean('thermal');
cmap = colormap("jet");

bypass_edge_together = together_plot(:, together_plot(5, :)==0); 
bypass_tip_together = together_plot(:, together_plot(5, :)==1);
pole_vaulting_together = together_plot(:, together_plot(5, :)==2); 
trapped_together = together_plot(:, together_plot(5, :)==3);

scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond'); hold on  % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'square'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', '^'); hold on % for legend only

scatter(trapped_together(6, :), trapped_together(4, :), 220, trapped_together(3, :), 'Filled', 'diamond','MarkerEdgeColor','k'); hold on 
scatter(bypass_edge_together(6, :), bypass_edge_together(4, :), 200, bypass_edge_together(3, :), 'Filled','MarkerEdgeColor','k'); hold on
scatter(bypass_tip_together(6, :), bypass_tip_together(4, :), 200, bypass_tip_together(3, :), 'Filled', 'square','MarkerEdgeColor','k'); hold on
scatter(pole_vaulting_together(6, :), pole_vaulting_together(4, :), 220, pole_vaulting_together(3, :), 'Filled', '^','MarkerEdgeColor','k'); hold on

% cmap(size(together_plot,2)); 
hcb=colorbar; caxis([0 1])
title(hcb,'$Initial\ position\ (y_0/h_{\rm obs})$','FontSize', 16,'Interpreter', 'latex'); grid on
set(gca,'FontSize',16);
% xlabel('$max\left|U_0-U(t)\right|/U_0$','FontSize', 22,'Interpreter', 'latex');
% xlabel('$Contact\ probability\ (disturbed\ layer = 40\mu{m}) $','FontSize', 22,'Interpreter', 'latex');  % for interaction3
% xlabel('$Normalized\ direct\ contact\ duration$','FontSize', 22,'Interpreter', 'latex');  
xlabel('$Normalized\ perturbed\ duration$','FontSize', 22,'Interpreter', 'latex'); 
% xlabel('$Interaction\ index$','FontSize', 22,'Interpreter', 'latex'); 
ylabel('$Deviation\ (\delta/h_{\rm obs})$','FontSize', 22,'Interpreter', 'latex');
legend({'Trapping','Below','Above','Pole-vaulting'}, 'Location', 'northeast','FontSize', 14,'Interpreter', 'latex')
% title('$0.5<L<1$','FontSize', 22,'Interpreter', 'latex')
% xlim([0 12]); %ylim([-0.4 0.8])
% f=gcf;
% exportgraphics(f,'F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\Figures\about interaction index\perturbed duration (0o30).png','Resolution',100)







%% # Contact infromation vs initial condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;
xlsfile = readcell(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\results_2023_02_23_Zhibo_DetectContactbyForce.xlsx'],'Sheet','Sheet1' ...
    ,'NumHeaderLines',1);
mask = cellfun(@ismissing, xlsfile); xlsfile(mask) = {nan};
together_plot = [cell2mat(xlsfile(:, 1:3)), cell2mat(xlsfile(:, 5)), ...
    cell2mat(xlsfile(:, 9:10)), cell2mat(xlsfile(:, 4)), cell2mat(xlsfile(:, 12)), ...
    cell2mat(xlsfile(:, 16))]; 
together_plot = together_plot';
prompt = {'The lower bound of the initial angle:', 'The upper bound of the initial angle:', ...
    'The lower bound of the contour length:','The upper bound of the contour length:'...
    'The lower bound of the initial position:','The upper bound of the initial position:'};
definput = {'-10', '10', 'nan', 'nan', 'nan', 'nan'};
answer = inputdlg(prompt, 'Input (please input NaN if there is no bound)', [1 35] , definput);

%%% assign the values 
range_chi0_low = str2double(answer{1,1}); if isnan(range_chi0_low); range_chi0_low = -91; end
range_chi0_up = str2double(answer{2,1}); if isnan(range_chi0_up); range_chi0_up = 91; end
range_L_low = str2double(answer{3,1}); if isnan(range_L_low); range_L_low = 0; end
range_L_up = str2double(answer{4,1}); if isnan(range_L_up); range_L_up = 10; end
range_y0_low = str2double(answer{5,1}); if isnan(range_y0_low); range_y0_low = -10; end
range_y0_up = str2double(answer{6,1}); if isnan(range_y0_up); range_y0_up = 10; end
% chi_0
together_plot(:, together_plot(1, :) < range_chi0_low) = []; 
together_plot(:, together_plot(1, :) > range_chi0_up) = [];
% L
together_plot(:, together_plot(2, :) < range_L_low) = [];
together_plot(:, together_plot(2, :) > range_L_up) = [];
% y_0
together_plot(:, together_plot(3, :) < range_y0_low) = [];  
together_plot(:, together_plot(3, :) > range_y0_up) = [];

together_plot([5:6, 8], or(isnan(together_plot(9, :)), together_plot(9, :) == 0)) = nan;  % use this for vector map plotting.
% together_plot(:, or(isnan(together_plot(9, :)), together_plot(9, :) == 0)) = []; % use this for the other plottings.
% choose the contact cases where the fiber ends touch the obstacle edge and set others to be NaN.

bypass_edge_together = together_plot(:, together_plot(4, :)==0); 
bypass_tip_together = together_plot(:, together_plot(4, :)==1);
pole_vaulting_together = together_plot(:, together_plot(4, :)==2); 
trapped_together = together_plot(:, together_plot(4, :)==3);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot theta_c vs y_c (with dynamics) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('color', 'w'); set(gcf, 'Position', [100 100 420 315]);
% for legend 
plot(nan, nan, 'o','MarkerSize', 8,'MarkerEdgeColor','k','MarkerFaceColor','red');
plot(nan, nan,  'square','MarkerSize', 8,'MarkerEdgeColor','k','MarkerFaceColor', [0 .5 0]); 
plot(nan, nan,  '^','MarkerSize', 8,'MarkerEdgeColor','k','MarkerFaceColor','blue');  
plot(nan, nan, 'diamond','MarkerSize', 8,'MarkerEdgeColor','k','MarkerFaceColor',[0.92, 0.70, 0.22]);
% for plot
plot(bypass_edge_together(6, :), bypass_edge_together(5, :), 'o','MarkerSize', 18,'MarkerEdgeColor','k','MarkerFaceColor','red'); hold on
plot(bypass_tip_together(6, :), bypass_tip_together(5, :),  'square','MarkerSize', 20,'MarkerEdgeColor','k','MarkerFaceColor', [0 .5 0]); hold on
plot(pole_vaulting_together(6, :), pole_vaulting_together(5, :),  '^','MarkerSize', 18,'MarkerEdgeColor','k','MarkerFaceColor','blue'); hold on
plot(trapped_together(6, :), trapped_together(5, :), 'diamond','MarkerSize', 18,'MarkerEdgeColor','k','MarkerFaceColor',[0.92, 0.70, 0.22]); hold on 

xlabel('$\theta_{\rm c}$','FontSize', 14,'Interpreter', 'latex'); 
ylabel('$y_{\rm c}/h_{\rm obs}$','FontSize', 14,'Interpreter', 'latex');
% title_txt = ['$-90 < \theta_0 < 90$'];
% title(title_txt,'FontSize', 18,'Interpreter', 'latex');
xlim([-50 50]); ylim([0 1]);
xticks([-45 -30 -15 0 15 30 45])
set(gca,'TickLabelInterpreter','latex','Fontsize',14);
grid on
% legend({'Trapping','Below','Above','Pole-vaulting'}, 'Location', 'northwest','FontSize', 14,'Interpreter', 'latex')

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\Simu data 2023-02-23\thetaC-yC_dyn_simu_from_far_theta_0m10to10_L1.png'],'Resolution',100)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot theta_c vs y_c (without dynamics) ALL DATA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('color', 'w'); set(gcf, 'Position', [100 100 520 390]);

% for plot
plot(together_plot(6, :), together_plot(5, :), 'o','MarkerSize', 10,'MarkerEdgeColor','k','MarkerFaceColor','m'); hold on 

xlabel('$\theta_{\rm c}$','FontSize', 18,'Interpreter', 'latex'); 
ylabel('$y_{\rm c}/h_{\rm obs}$','FontSize', 18,'Interpreter', 'latex');
xlim([-90 90]); ylim([0 1])
set(gca,'TickLabelInterpreter','latex','Fontsize',18);

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\Simu data 2023-02-23\thetaC-yC_simu_from_far_theta_0m10to10_alldata.png'],'Resolution',100)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot theta_0 vs y_0 (without dynamics) which contact the obstacle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('color', 'w'); set(gcf, 'Position', [100 100 520 390]);

% choose the contact cases
together_plot_contact = together_plot;
together_plot_contact(:, isnan(together_plot_contact(5, :))) = [];

% for plot
plot(together_plot_contact(1, :), together_plot_contact(3, :), ...
    'o','MarkerSize', 10,'MarkerEdgeColor','k','MarkerFaceColor','m'); hold on 

xlabel('$\theta_0$','FontSize', 18,'Interpreter', 'latex'); 
ylabel('$y_0/h_{\rm obs}$','FontSize', 18,'Interpreter', 'latex');
xlim([-10 10]); ylim([0 1])
set(gca,'TickLabelInterpreter','latex','Fontsize',18);

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\Simu data 2023-02-23\theta0-y0_from_far_contactcases.png'],'Resolution',100)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot theta_0 vs y_0 (with dynamics) which contact the obstacle 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('color', 'w'); set(gcf, 'Position', [100 100 600 450]);
% for legend
plot(nan, nan, 'diamond','MarkerSize', 8,'MarkerEdgeColor','k','MarkerFaceColor',[0.92, 0.70, 0.22]); 
plot(nan, nan, 'o','MarkerSize', 8,'MarkerEdgeColor','k','MarkerFaceColor','red');
plot(nan, nan,  'square','MarkerSize', 8,'MarkerEdgeColor','k','MarkerFaceColor', [0 .5 0]); 
plot(nan, nan,  '^','MarkerSize', 8,'MarkerEdgeColor','k','MarkerFaceColor','blue'); 
% choose the contact cases
together_plot_contact = together_plot;
together_plot_contact(:, isnan(together_plot_contact(5, :))) = [];
bypass_edge_together_contact = together_plot_contact(:, together_plot_contact(4, :)==0); 
bypass_tip_together_contact = together_plot_contact(:, together_plot_contact(4, :)==1);
pole_vaulting_together_contact = together_plot_contact(:, together_plot_contact(4, :)==2); 
trapped_together_contact = together_plot_contact(:, together_plot_contact(4, :)==3);

% for plot
plot(trapped_together_contact(1, :), trapped_together_contact(3, :), ...
    'diamond','MarkerSize', 10,'MarkerEdgeColor','k','MarkerFaceColor',[0.92, 0.70, 0.22]); hold on 
plot(bypass_edge_together_contact(1, :), bypass_edge_together_contact(3, :), ...
    'o','MarkerSize', 10,'MarkerEdgeColor','k','MarkerFaceColor','red'); hold on
plot(bypass_tip_together_contact(1, :), bypass_tip_together_contact(3, :), ...
    'square','MarkerSize', 10,'MarkerEdgeColor','k','MarkerFaceColor', [0 .5 0]); hold on
plot(pole_vaulting_together_contact(1, :), pole_vaulting_together_contact(3, :), ...
    '^','MarkerSize', 10,'MarkerEdgeColor','k','MarkerFaceColor','blue'); hold on

xlabel('$\theta_0$','FontSize', 18,'Interpreter', 'latex'); 
ylabel('$y_0/h_{\rm obs}$','FontSize', 18,'Interpreter', 'latex');
% title_txt = ['$-90 < \theta_0 < 90$'];
% title(title_txt,'FontSize', 18,'Interpreter', 'latex');
xlim([-10 10]); ylim([0 1])
set(gca,'TickLabelInterpreter','latex','Fontsize',18);
% xticks([-45 -30 -15 0 15 30 45])
% legend({'Trapping','Below','Above','Pole-vaulting'}, 'Location', 'northwest','FontSize', 14,'Interpreter', 'latex')

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\Simu data 2023-02-23\theta0-y0_dyn_simu_from_far_contactcases.png'],'Resolution',100)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot theta_c vs theta_0 (with color y_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('color', 'w'); set(gcf, 'Position', [100 100 760 570]);

% together_plot_less_y0 = together_plot(:, or(or(or(or(or(or(together_plot(3, :)==0.125, together_plot(3, :)==0.250), ...
% together_plot(3, :)==0.375), together_plot(3, :)==0.5), together_plot(3, :)==0.625), ...
% together_plot(3, :)==0.75),together_plot(3, :)==0.875)); % choose y_0

% together_plot_less_y0 = together_plot(:, or(or(or(together_plot(3, :)==0.125, ...
% together_plot(3, :)==0.375), together_plot(3, :)==0.625), together_plot(3, :)==0.875));  % choose y_0

together_plot_less_y0 = together_plot;  

L1_together = together_plot_less_y0(:, together_plot_less_y0(2, :)==0.5); % classify contour length and indicate by symbols
L2_together = together_plot_less_y0(:, together_plot_less_y0(2, :)==0.6); 
L3_together = together_plot_less_y0(:, together_plot_less_y0(2, :)==0.7); 
L4_together = together_plot_less_y0(:, together_plot_less_y0(2, :)==0.8); 
L5_together = together_plot_less_y0(:, together_plot_less_y0(2, :)==0.9); 
L6_together = together_plot_less_y0(:, together_plot_less_y0(2, :)==1); 
L7_together = together_plot_less_y0(:, together_plot_less_y0(2, :)==1.2);
L8_together = together_plot_less_y0(:, together_plot_less_y0(2, :)==1.4);

% for legend
scatter(nan, nan, 1, nan, 'Filled', 'diamond','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on  % for legend only
scatter(nan, nan, 1, nan, 'Filled', 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only
scatter(nan, nan, 1, nan, 'Filled', 'square','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only
scatter(nan, nan, 1, nan, 'Filled', '^','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only

scatter(L1_together(1, :), L1_together(6, :), 100, L1_together(3, :), 'Filled', 'diamond','MarkerEdgeColor','k'); hold on 
scatter(L3_together(1, :), L3_together(6, :), 100, L3_together(3, :), 'Filled','o', 'MarkerEdgeColor','k'); hold on
scatter(L5_together(1, :), L5_together(6, :), 100, L5_together(3, :), 'Filled', 'square','MarkerEdgeColor','k'); hold on
scatter(L7_together(1, :), L7_together(6, :), 100, L7_together(3, :), 'Filled', '^','MarkerEdgeColor','k'); hold on
% scatter(L2_together(1, :), L2_together(6, :), 100, L2_together(3, :), 'diamond','LineWidth',2); hold on 
% scatter(L4_together(1, :), L4_together(6, :), 100, L4_together(3, :), 'o', 'LineWidth',2); hold on
% scatter(L6_together(1, :), L6_together(6, :), 100, L6_together(3, :), 'square','LineWidth',2); hold on
% scatter(L8_together(1, :), L8_together(6, :), 100, L8_together(3, :), '^','LineWidth',2); hold on

% cmap(size(together_plot,2)); 
hcb=colorbar('Ticks', 0.1:0.1:0.6); caxis([0.1 0.6]); colormap jet
hcb.FontSize = 16;
title(hcb,'$y_0/h_{\rm obs}$','FontSize', 24,'Interpreter', 'latex'); grid on
xlabel('$\theta_0$','FontSize', 24,'Interpreter', 'latex'); 
ylabel('$\theta_{\rm c}$','FontSize', 24,'Interpreter', 'latex');

% title_txt = ['$-10 < \theta_0 < 10$'];
% title(title_txt,'FontSize', 18,'Interpreter', 'latex');
% xlim([-90 90]); ylim([-0.1 1.1]);
legend({'$L/l_{\rm obs}=0.5$','$L/l_{\rm obs}=0.7$','$L/l_{\rm obs}=0.9$','$L/l_{\rm obs}=1.2$'}, 'Location', 'northwest','FontSize', 16,'Interpreter', 'latex')

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\Simu data 2023-02-23\thetac-theta0_color-y0.png'],'Resolution',100)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot theta_c vs theta_0 (with color y_c) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('color', 'w'); set(gcf, 'Position', [100 100 760 570]);

L1_together = together_plot(:, together_plot(2, :)==0.5); % classify contour length and indicate by symbols
L2_together = together_plot(:, together_plot(2, :)==0.6); 
L3_together = together_plot(:, together_plot(2, :)==0.7); 
L4_together = together_plot(:, together_plot(2, :)==0.8); 
L5_together = together_plot(:, together_plot(2, :)==0.9); 
L6_together = together_plot(:, together_plot(2, :)==1); 
L7_together = together_plot(:, together_plot(2, :)==1.2);
L8_together = together_plot(:, together_plot(2, :)==1.4);

% for legend
scatter(nan, nan, 1, nan, 'Filled', 'diamond','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on  % for legend only
scatter(nan, nan, 1, nan, 'Filled', 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only
scatter(nan, nan, 1, nan, 'Filled', 'square','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only
scatter(nan, nan, 1, nan, 'Filled', '^','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only

scatter(L1_together(1, :), L1_together(6, :), 100, L1_together(5, :), 'Filled', 'diamond','MarkerEdgeColor','k'); hold on 
scatter(L3_together(1, :), L3_together(6, :), 100, L3_together(5, :), 'Filled','o', 'MarkerEdgeColor','k'); hold on
scatter(L5_together(1, :), L5_together(6, :), 100, L5_together(5, :), 'Filled', 'square','MarkerEdgeColor','k'); hold on
scatter(L7_together(1, :), L7_together(6, :), 100, L7_together(5, :), 'Filled', '^','MarkerEdgeColor','k'); hold on
% scatter(L2_together(1, :), L2_together(6, :), 100, L2_together(5, :), 'diamond','LineWidth',2); hold on 
% scatter(L4_together(1, :), L4_together(6, :), 100, L4_together(5, :), 'o', 'LineWidth',2); hold on
% scatter(L6_together(1, :), L6_together(6, :), 100, L6_together(5, :), 'square','LineWidth',2); hold on
% scatter(L8_together(1, :), L8_together(6, :), 100, L8_together(5, :), '^','LineWidth',2); hold on

% cmap(size(together_plot,2)); 
hcb=colorbar('Ticks', -0.1:0.2:1.1); caxis([-0.1 1.1]); colormap jet
hcb.FontSize = 16;
title(hcb,'$y_{\rm c}/h_{\rm obs}$','FontSize', 24,'Interpreter', 'latex'); grid on
xlabel('$\theta_0$','FontSize', 24,'Interpreter', 'latex'); 
ylabel('$\theta_{\rm c}$','FontSize', 24,'Interpreter', 'latex');

% title_txt = ['$-10 < \theta_0 < 10$'];
% title(title_txt,'FontSize', 18,'Interpreter', 'latex');
% xlim([-90 90]); ylim([-0.1 1.1]);
legend({'$L/l_{\rm obs}=0.5$','$L/l_{\rm obs}=0.7$','$L/l_{\rm obs}=0.9$','$L/l_{\rm obs}=1.2$'}, 'Location', 'northwest','FontSize', 16,'Interpreter', 'latex')

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\Simu data 2023-02-23\thetac-theta0_color-yc.png'],'Resolution',100)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot (theta_c-theta_0) vs theta_0 (with color y_0) 
% (theta_c-theta_0): the rotation angle before the first contact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('color', 'w'); set(gcf, 'Position', [100 100 760 570]);

% together_plot_less_y0 = together_plot(:, or(or(or(or(or(or(together_plot(3, :)==0.125, together_plot(3, :)==0.250), ...
% together_plot(3, :)==0.375), together_plot(3, :)==0.5), together_plot(3, :)==0.625), ...
% together_plot(3, :)==0.75),together_plot(3, :)==0.875)); % choose y_0

% together_plot_less_y0 = together_plot(:, or(or(or(together_plot(3, :)==0.125, ...
% together_plot(3, :)==0.375), together_plot(3, :)==0.625), together_plot(3, :)==0.875));  % choose y_0

together_plot_less_y0 = together_plot;  

L1_together = together_plot_less_y0(:, together_plot_less_y0(2, :)==0.5); % classify contour length and indicate by symbols
L2_together = together_plot_less_y0(:, together_plot_less_y0(2, :)==0.6); 
L3_together = together_plot_less_y0(:, together_plot_less_y0(2, :)==0.7); 
L4_together = together_plot_less_y0(:, together_plot_less_y0(2, :)==0.8); 
L5_together = together_plot_less_y0(:, together_plot_less_y0(2, :)==0.9); 
L6_together = together_plot_less_y0(:, together_plot_less_y0(2, :)==1); 
L7_together = together_plot_less_y0(:, together_plot_less_y0(2, :)==1.2);
L8_together = together_plot_less_y0(:, together_plot_less_y0(2, :)==1.4);

% for legend
scatter(nan, nan, 1, nan, 'Filled', 'diamond','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on  % for legend only
scatter(nan, nan, 1, nan, 'Filled', 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only
scatter(nan, nan, 1, nan, 'Filled', 'square','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only
scatter(nan, nan, 1, nan, 'Filled', '^','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only

scatter(L1_together(1, :), L1_together(8, :), 100, L1_together(3, :), 'Filled', 'diamond','MarkerEdgeColor','k'); hold on 
scatter(L3_together(1, :), L3_together(8, :), 100, L3_together(3, :), 'Filled','o', 'MarkerEdgeColor','k'); hold on
scatter(L5_together(1, :), L5_together(8, :), 100, L5_together(3, :), 'Filled', 'square','MarkerEdgeColor','k'); hold on
scatter(L7_together(1, :), L7_together(8, :), 100, L7_together(3, :), 'Filled', '^','MarkerEdgeColor','k'); hold on
% scatter(L2_together(1, :), L2_together(8, :), 100, L2_together(3, :), 'diamond','LineWidth',2); hold on 
% scatter(L4_together(1, :), L4_together(8, :), 100, L4_together(3, :), 'o', 'LineWidth',2); hold on
% scatter(L6_together(1, :), L6_together(8, :), 100, L6_together(3, :), 'square','LineWidth',2); hold on
% scatter(L8_together(1, :), L8_together(8, :), 100, L8_together(3, :), '^','LineWidth',2); hold on

% cmap(size(together_plot,2)); 
hcb=colorbar('Ticks', 0.1:0.1:0.6); caxis([0.1 0.6]); colormap jet
hcb.FontSize = 16;
title(hcb,'$y_0/h_{\rm obs}$','FontSize', 24,'Interpreter', 'latex'); grid on
xlabel('$\theta_0$','FontSize', 24,'Interpreter', 'latex'); 
ylabel('$\theta_{\rm c}-\theta_0$','FontSize', 24,'Interpreter', 'latex');

% title_txt = ['$-10 < \theta_0 < 10$'];
% title(title_txt,'FontSize', 18,'Interpreter', 'latex');
% xlim([-90 90]); ylim([-0.1 1.1]);
legend({'$L/l_{\rm obs}=0.5$','$L/l_{\rm obs}=0.7$','$L/l_{\rm obs}=0.9$','$L/l_{\rm obs}=1.2$'}, 'Location', 'northwest','FontSize', 16,'Interpreter', 'latex')

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\Simu data 2023-02-23\deltatheta-theta0_color-y0.png'],'Resolution',100)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot y_c vs y_0 (with color theta_0) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('color', 'w'); set(gcf, 'Position', [100 100 760 570]);

% together_plot_less_theta0 = together_plot(:, or(or(or(or(or(together_plot(1, :)==-60, together_plot(1, :)==-30), ...
% together_plot(1, :)==0), together_plot(1, :)==30), together_plot(1, :)==60), together_plot(1, :)==90)); % choose theta_0

% together_plot_less_theta0 = together_plot(:, or(or(or(or(or(or(or(or(together_plot(1, :)==-10, together_plot(1, :)==-7.5), ...
% together_plot(1, :)==-5), together_plot(1, :)==-2.5), together_plot(1, :)==0), ...
% together_plot(1, :)==2.5), together_plot(1, :)==5), together_plot(1, :)==7.5), together_plot(1, :)==10)); % choose theta_0

together_plot_less_theta0 = together_plot;  

L1_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==0.5); % classify contour length and indicate by symbols
L2_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==0.6); 
L3_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==0.7); 
L4_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==0.8); 
L5_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==0.9); 
L6_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==1); 
L7_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==1.2);
L8_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==1.4);

% for legend
scatter(nan, nan, 1, nan, 'Filled', 'diamond','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on  % for legend only
scatter(nan, nan, 1, nan, 'Filled', 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only
scatter(nan, nan, 1, nan, 'Filled', 'square','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only
scatter(nan, nan, 1, nan, 'Filled', '^','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only

scatter(L1_together(3, :), L1_together(5, :), 100, L1_together(1, :), 'Filled', 'diamond','MarkerEdgeColor','k'); hold on 
scatter(L3_together(3, :), L3_together(5, :), 100, L3_together(1, :), 'Filled','o', 'MarkerEdgeColor','k'); hold on
scatter(L5_together(3, :), L5_together(5, :), 100, L5_together(1, :), 'Filled', 'square','MarkerEdgeColor','k'); hold on
scatter(L7_together(3, :), L7_together(5, :), 100, L7_together(1, :), 'Filled', '^','MarkerEdgeColor','k'); hold on
% scatter(L2_together(3, :), L2_together(5, :), 100, L2_together(1, :), 'diamond','LineWidth',2); hold on 
% scatter(L4_together(3, :), L4_together(5, :), 100, L4_together(1, :), 'o', 'LineWidth',2); hold on
% scatter(L6_together(3, :), L6_together(5, :), 100, L6_together(1, :), 'square','LineWidth',2); hold on
% scatter(L8_together(3, :), L8_together(5, :), 100, L8_together(1, :), '^','LineWidth',2); hold on

hcb=colorbar('Ticks', -10:2.5:10); caxis([-10 10]); colormap jet
hcb.FontSize = 16;
title(hcb,'$\theta_0$','FontSize', 24,'Interpreter', 'latex'); grid on

xlabel('$y_0/h_{\rm obs}$','FontSize', 24,'Interpreter', 'latex'); 
ylabel('$y_{\rm c}/h_{\rm obs}$','FontSize', 24,'Interpreter', 'latex');
% title_txt = ['$-10 < \theta_0 < 10$'];
% title(title_txt,'FontSize', 18,'Interpreter', 'latex');
ylim([-0.1 1.1]);
legend({'$L/l_{\rm obs}=0.5$','$L/l_{\rm obs}=0.7$','$L/l_{\rm obs}=0.9$','$L/l_{\rm obs}=1.2$'}, 'Location', 'northwest','FontSize', 16,'Interpreter', 'latex')

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\Simu data 2023-02-23\yc-y0_color-theta_0.png'],'Resolution',100)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot theta_0, y_0, L and colorcoded (theta_c-theta_0)
% (theta_c-theta_0): the rotation angle before the first contact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('color', 'w'); set(gcf, 'Position', [100 100 1000 600]);

[XX, YY, ZZ] = meshgrid(unique(together_plot(1, :)), unique(together_plot(3, :)), ...
    unique(together_plot(2, :)));
scatter3(XX(:), YY(:), ZZ(:), 100, ones(size(ZZ(:), 1), 1), 'o', ...
    'MarkerEdgeColor', [.8 .8 .8]); hold on 

scatter3(together_plot(1, :), together_plot(3, :), together_plot(2, :), 100, ... 
together_plot(8, :), 'Filled', 'o','MarkerEdgeColor','k');

% [XX, YY, ZZ] = meshgrid(unique(together_plot(1, :)), unique(together_plot(3, :)), unique(together_plot(2, :)));
% VV = double(logical(XX));
% hh1 = slice(XX, YY, ZZ, VV, [], 0.2, []);
% hh2 = slice(XX, YY, ZZ, VV, [], 0.4, []);
% 
% set(hh1, 'EdgeColor','none','FaceColor','k','FaceAlpha', 0.2);
% set(hh2, 'EdgeColor','none','FaceColor','k','FaceAlpha', 0.4);

hcb=colorbar('Ticks', -40:20:40); caxis([-40 40]); cmocean('balance')
hcb.FontSize = 14;
ylabel(hcb,'$Fiber\ rotation\ angle\ before\ contact$','FontSize', 18,'Interpreter', 'latex'); grid on
xlabel('$\theta_0$','FontSize', 22,'Interpreter', 'latex'); 
ylabel('$y_0/h_{\rm obs}$','FontSize', 22,'Interpreter', 'latex'); yticks([0.2 0.4])
zlabel('$L/l_{\rm obs}$','FontSize', 22,'Interpreter', 'latex');
set(gca,'TickLabelInterpreter','latex');

view(-12,5)

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\Simu data 2023-02-23\theta_0-y_0-L_vs_RotationAngle.png'],'Resolution',100)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot theta_0, y_0, L and colorcoded (y_c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('color', 'w'); set(gcf, 'Position', [100 100 1000 600]);

[XX, YY, ZZ] = meshgrid(unique(together_plot(1, :)), unique(together_plot(3, :)), ...
    unique(together_plot(2, :)));
scatter3(XX(:), YY(:), ZZ(:), 100, ones(size(ZZ(:), 1), 1), 'o', ...
    'MarkerEdgeColor', [.8 .8 .8]); hold on 

scatter3(together_plot(1, :), together_plot(3, :), together_plot(2, :), 100, ... 
together_plot(5, :), 'Filled', 'o','MarkerEdgeColor','k'); 

% [XX, YY, ZZ] = meshgrid(unique(together_plot(1, :)), unique(together_plot(3, :)), unique(together_plot(2, :)));
% VV = double(logical(XX));
% hh1 = slice(XX, YY, ZZ, VV, [], 0.2, []);
% hh2 = slice(XX, YY, ZZ, VV, [], 0.4, []);
% 
% set(hh1, 'EdgeColor','none','FaceColor','k','FaceAlpha', 0.2);
% set(hh2, 'EdgeColor','none','FaceColor','k','FaceAlpha', 0.4);

hcb=colorbar('Ticks', 0:0.2:1); caxis([0 1]); cmocean('balance')
hcb.FontSize = 14;
ylabel(hcb,'$y_{\rm c}/h_{\rm obs}$','FontSize', 18,'Interpreter', 'latex'); grid on
xlabel('$\theta_0$','FontSize', 22,'Interpreter', 'latex'); 
ylabel('$y_0/h_{\rm obs}$','FontSize', 22,'Interpreter', 'latex'); yticks([0.2 0.4])
zlabel('$L/l_{\rm obs}$','FontSize', 22,'Interpreter', 'latex');
set(gca,'TickLabelInterpreter','latex');

view(-12,5)

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\Simu data 2023-02-23\theta_0-y_0-L_vs_y_c.png'],'Resolution',100)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot y_c_CoM vs y_0 (with color theta_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('color', 'w'); set(gcf, 'Position', [100 100 760 570]);

together_plot_less_theta0 = together_plot;  

L1_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==0.5); % classify contour length and indicate by symbols
L2_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==0.6); 
L3_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==0.7); 
L4_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==0.8); 
L5_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==0.9); 
L6_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==1); 
L7_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==1.2);
L8_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==1.4);

% for legend
scatter(nan, nan, 1, nan, 'Filled', 'diamond','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on  % for legend only
scatter(nan, nan, 1, nan, 'Filled', 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only
scatter(nan, nan, 1, nan, 'Filled', 'square','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only
scatter(nan, nan, 1, nan, 'Filled', '^','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only

scatter(L1_together(3, :), L1_together(5, :)-0.5/2*sind(L1_together(6, :)), 100, L1_together(1, :), 'Filled', 'diamond','MarkerEdgeColor','k'); hold on 
scatter(L3_together(3, :), L3_together(5, :)-0.7/2*sind(L3_together(6, :)), 100, L3_together(1, :), 'Filled','o', 'MarkerEdgeColor','k'); hold on
scatter(L5_together(3, :), L5_together(5, :)-0.9/2*sind(L5_together(6, :)), 100, L5_together(1, :), 'Filled', 'square','MarkerEdgeColor','k'); hold on
scatter(L7_together(3, :), L7_together(5, :)-1.2/2*sind(L7_together(6, :)), 100, L7_together(1, :), 'Filled', '^','MarkerEdgeColor','k'); hold on
% scatter(L2_together(3, :), L2_together(5, :)-0.6/2*sind(L2_together(6, :)), 100, L2_together(1, :), 'diamond','LineWidth',2); hold on 
% scatter(L4_together(3, :), L4_together(5, :)-0.8/2*sind(L4_together(6, :)), 100, L4_together(1, :), 'o', 'LineWidth',2); hold on
% scatter(L6_together(3, :), L6_together(5, :)-1.0/2*sind(L6_together(6, :)), 100, L6_together(1, :), 'square','LineWidth',2); hold on
% scatter(L8_together(3, :), L8_together(5, :)-1.4/2*sind(L8_together(6, :)), 100, L8_together(1, :), '^','LineWidth',2); hold on

hcb=colorbar('Ticks', -10:2.5:10); caxis([-10 10]); colormap jet
hcb.FontSize = 16;
title(hcb,'$\theta_0$','FontSize', 24,'Interpreter', 'latex'); grid on

xlabel('$y_0/h_{\rm obs}$','FontSize', 24,'Interpreter', 'latex'); 
ylabel('${y_{\rm c}}^{CoM}/h_{\rm obs}$','FontSize', 24,'Interpreter', 'latex');
% title_txt = ['$-10 < \theta_0 < 10$'];
% title(title_txt,'FontSize', 18,'Interpreter', 'latex');
ylim([-0.1 1.1]);
legend({'$L/l_{\rm obs}=0.5$','$L/l_{\rm obs}=0.7$','$L/l_{\rm obs}=0.9$','$L/l_{\rm obs}=1.2$'}, 'Location', 'northwest','FontSize', 16,'Interpreter', 'latex')

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\Simu data 2023-02-23\yc_COM-y0_color-theta_0.png'],'Resolution',100)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot y_{\rm c} vs y_0_HEAD (with color theta_0) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('color', 'w'); set(gcf, 'Position', [100 100 760 570]);

% together_plot_less_theta0 = together_plot(:, or(or(or(or(or(together_plot(1, :)==-60, together_plot(1, :)==-30), ...
% together_plot(1, :)==0), together_plot(1, :)==30), together_plot(1, :)==60), together_plot(1, :)==90)); % choose theta_0

% together_plot_less_theta0 = together_plot(:, or(or(or(or(or(or(or(or(together_plot(1, :)==-10, together_plot(1, :)==-7.5), ...
% together_plot(1, :)==-5), together_plot(1, :)==-2.5), together_plot(1, :)==0), ...
% together_plot(1, :)==2.5), together_plot(1, :)==5), together_plot(1, :)==7.5), together_plot(1, :)==10)); % choose theta_0

together_plot_less_theta0 = together_plot;  

L1_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==0.5); % classify contour length and indicate by symbols
L2_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==0.6); 
L3_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==0.7); 
L4_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==0.8); 
L5_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==0.9); 
L6_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==1); 
L7_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==1.2);
L8_together = together_plot_less_theta0(:, together_plot_less_theta0(2, :)==1.4);

% for legend
scatter(nan, nan, 1, nan, 'Filled', 'diamond','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on  % for legend only
scatter(nan, nan, 1, nan, 'Filled', 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only
scatter(nan, nan, 1, nan, 'Filled', 'square','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only
scatter(nan, nan, 1, nan, 'Filled', '^','MarkerEdgeColor','k', 'MarkerFaceColor',[.7 .7 .7]); hold on % for legend only

scatter(L1_together(3, :)+0.5/2*sind(L1_together(1, :)), L1_together(5, :), 100, L1_together(1, :), 'Filled', 'diamond','MarkerEdgeColor','k'); hold on 
scatter(L3_together(3, :)+0.7/2*sind(L3_together(1, :)), L3_together(5, :), 100, L3_together(1, :), 'Filled','o', 'MarkerEdgeColor','k'); hold on
scatter(L5_together(3, :)+0.9/2*sind(L5_together(1, :)), L5_together(5, :), 100, L5_together(1, :), 'Filled', 'square','MarkerEdgeColor','k'); hold on
scatter(L7_together(3, :)+1.2/2*sind(L7_together(1, :)), L7_together(5, :), 100, L7_together(1, :), 'Filled', '^','MarkerEdgeColor','k'); hold on
% scatter(L2_together(3, :)+0.6/2*sind(L2_together(1, :)), L2_together(5, :), 100, L2_together(1, :), 'diamond','LineWidth',2); hold on 
% scatter(L4_together(3, :)+0.8/2*sind(L4_together(1, :)), L4_together(5, :), 100, L4_together(1, :), 'o', 'LineWidth',2); hold on
% scatter(L6_together(3, :)+1.0/2*sind(L6_together(1, :)), L6_together(5, :), 100, L6_together(1, :), 'square','LineWidth',2); hold on
% scatter(L8_together(3, :)+1.4/2*sind(L8_together(1, :)), L8_together(5, :), 100, L8_together(1, :), '^','LineWidth',2); hold on

hcb=colorbar('Ticks', -10:2.5:10); caxis([-10 10]); colormap jet
hcb.FontSize = 18; hcb.TickLabelInterpreter = 'latex';
title(hcb,'$\theta_0$','FontSize', 18,'Interpreter', 'latex'); grid on

xlabel('$y_0^{HEAD}/h_{\rm obs}$','FontSize', 18,'Interpreter', 'latex'); 
ylabel('$y_{\rm c}/h_{\rm obs}$','FontSize', 18,'Interpreter', 'latex');
% title_txt = ['$-10 < \theta_0 < 10$'];
% title(title_txt,'FontSize', 18,'Interpreter', 'latex');
ylim([0 1]);
set(gca,'TickLabelInterpreter','latex','Fontsize',18);
legend({'$L/l_{\rm obs}=0.5$','$L/l_{\rm obs}=0.7$','$L/l_{\rm obs}=0.9$','$L/l_{\rm obs}=1.2$'}, 'Location', 'northwest','FontSize', 16,'Interpreter', 'latex')

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\Simu data 2023-02-23\yc-y0_HEAD_color-theta_0.png'],'Resolution',100)





% %%%%%%%%%%%% plot theta_c vs y_c (with initial position y_0) %%%%%%%%%%%%%%
% 
% figure('color', 'w'); set(gcf, 'Position', [100 100 1500 300]);
% % for legend
% scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond'); hold on  % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', 'square'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', '^'); hold on % for legend only
% 
% scatter(trapped_together(6, :), trapped_together(5, :), 150, trapped_together(3, :), 'Filled', 'diamond','MarkerEdgeColor','k'); hold on 
% scatter(bypass_edge_together(6, :), bypass_edge_together(5, :), 130, bypass_edge_together(3, :), 'Filled','MarkerEdgeColor','k'); hold on
% scatter(bypass_tip_together(6, :), bypass_tip_together(5, :), 130, bypass_tip_together(3, :), 'Filled', 'square','MarkerEdgeColor','k'); hold on
% scatter(pole_vaulting_together(6, :), pole_vaulting_together(5, :), 150, pole_vaulting_together(3, :), 'Filled', '^','MarkerEdgeColor','k'); hold on
% 
% % cmap(size(together_plot,2)); 
% hcb=colorbar; caxis([0 1]); colormap jet
% title(hcb,'$Initial\ position\ (y_0/h_{\rm obs})$','FontSize', 16,'Interpreter', 'latex'); grid on
% 
% xlabel('$\theta_c$','FontSize', 18,'Interpreter', 'latex'); 
% ylabel('$y_c$','FontSize', 18,'Interpreter', 'latex');
% title_txt = ['$-10 < \theta_0 < 10$'];
% title(title_txt,'FontSize', 18,'Interpreter', 'latex');
% xlim([-90 90]); ylim([-0.1 1.1]);
% legend({'Trapping','Below','Above','Pole-vaulting'}, 'Location', 'northwest','FontSize', 14,'Interpreter', 'latex')
% 
% % f=gcf;
% % exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
% %     '\Figures\about contact information vs initial condition\theta_0m10to10_theta_c-y_c_y_0.png'],'Resolution',100)
% 
% 
% 
% %%%%%%%%%%%% plot theta_c vs y_c (with initial position theta_0) %%%%%%%%%%%%%%
% 
% figure('color', 'w'); set(gcf, 'Position', [100 100 1500 300]);
% % for legend
% scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond'); hold on  % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', 'square'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', '^'); hold on % for legend only
% 
% scatter(trapped_together(6, :), trapped_together(5, :), 150, trapped_together(1, :), 'Filled', 'diamond','MarkerEdgeColor','k'); hold on 
% scatter(bypass_edge_together(6, :), bypass_edge_together(5, :), 130, bypass_edge_together(1, :), 'Filled','MarkerEdgeColor','k'); hold on
% scatter(bypass_tip_together(6, :), bypass_tip_together(5, :), 130, bypass_tip_together(1, :), 'Filled', 'square','MarkerEdgeColor','k'); hold on
% scatter(pole_vaulting_together(6, :), pole_vaulting_together(5, :), 150, pole_vaulting_together(1, :), 'Filled', '^','MarkerEdgeColor','k'); hold on
% 
% % cmap(size(together_plot,2)); 
% hcb=colorbar; caxis([-10 10]); colormap jet
% title(hcb,'$Initial\ angle\ ({\theta}_0)$','FontSize', 16,'Interpreter', 'latex'); grid on
% 
% xlabel('$\theta_c$','FontSize', 18,'Interpreter', 'latex'); 
% ylabel('$y_c$','FontSize', 18,'Interpreter', 'latex');
% title_txt = ['$-10 < \theta_0 < 10$'];
% title(title_txt,'FontSize', 18,'Interpreter', 'latex');
% xlim([-90 90]); ylim([-0.1 1.1]);
% legend({'Trapping','Below','Above','Pole-vaulting'}, 'Location', 'northwest','FontSize', 14,'Interpreter', 'latex')
% 
% % f=gcf;
% % exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
% %     '\Figures\about contact information vs initial condition\theta_0m10to10_theta_c-y_c_theta_0.png'],'Resolution',100)
% 
% 
% 
% %%%%%%%%%%%%%%% plot theta_c vs y_c (with contour length L) %%%%%%%%%%%%%%%%%
% 
% figure('color', 'w'); set(gcf, 'Position', [100 100 1500 300]);
% % for legend
% scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond'); hold on  % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', 'square'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', '^'); hold on % for legend only
% 
% scatter(trapped_together(6, :), trapped_together(5, :), 150, trapped_together(2, :), 'Filled', 'diamond','MarkerEdgeColor','k'); hold on 
% scatter(bypass_edge_together(6, :), bypass_edge_together(5, :), 130, bypass_edge_together(2, :), 'Filled','MarkerEdgeColor','k'); hold on
% scatter(bypass_tip_together(6, :), bypass_tip_together(5, :), 130, bypass_tip_together(2, :), 'Filled', 'square','MarkerEdgeColor','k'); hold on
% scatter(pole_vaulting_together(6, :), pole_vaulting_together(5, :), 150, pole_vaulting_together(2, :), 'Filled', '^','MarkerEdgeColor','k'); hold on
% 
% % cmap(size(together_plot,2)); 
% hcb=colorbar; caxis([0.5 1.5]); colormap jet
% title(hcb,'$Contour\ length\ (L/l_{\rm obs})$','FontSize', 16,'Interpreter', 'latex'); grid on
% 
% xlabel('$\theta_c$','FontSize', 18,'Interpreter', 'latex'); 
% ylabel('$y_c$','FontSize', 18,'Interpreter', 'latex');
% title_txt = ['$-10 < \theta_0 < 10$'];
% title(title_txt,'FontSize', 18,'Interpreter', 'latex');
% xlim([-90 90]); ylim([-0.1 1.1]);
% legend({'Trapping','Below','Above','Pole-vaulting'}, 'Location', 'northwest','FontSize', 14,'Interpreter', 'latex')
% 
% % f=gcf;
% % exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
% %     '\Figures\about contact information vs initial condition\theta_0m10to10_theta_c-y_c_L.png'],'Resolution',100)
% 
% 
% 
% %%%%%%%%%%%%%%% plot theta_c vs y_c (with deviation delta) %%%%%%%%%%%%%%%%%
% 
% figure('color', 'w'); set(gcf, 'Position', [100 100 1500 300]);
% % for legend
% scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond'); hold on  % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', 'square'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', '^'); hold on % for legend only
% 
% scatter(trapped_together(6, :), trapped_together(5, :), 150, trapped_together(7, :), 'Filled', 'diamond','MarkerEdgeColor','k'); hold on 
% scatter(bypass_edge_together(6, :), bypass_edge_together(5, :), 130, bypass_edge_together(7, :), 'Filled','MarkerEdgeColor','k'); hold on
% scatter(bypass_tip_together(6, :), bypass_tip_together(5, :), 130, bypass_tip_together(7, :), 'Filled', 'square','MarkerEdgeColor','k'); hold on
% scatter(pole_vaulting_together(6, :), pole_vaulting_together(5, :), 150, pole_vaulting_together(7, :), 'Filled', '^','MarkerEdgeColor','k'); hold on
% 
% % cmap(size(together_plot,2)); 
% hcb=colorbar; caxis([-0.6 0.6]); colormap jet
% title(hcb,'$Deviation\ (\delta)$','FontSize', 16,'Interpreter', 'latex'); grid on
% 
% xlabel('$\theta_c$','FontSize', 18,'Interpreter', 'latex'); 
% ylabel('$y_c$','FontSize', 18,'Interpreter', 'latex');
% title_txt = ['$-10 < \theta_0 < 10$'];
% title(title_txt,'FontSize', 18,'Interpreter', 'latex');
% xlim([-90 90]); ylim([-0.1 1.1]);
% legend({'Trapping','Below','Above','Pole-vaulting'}, 'Location', 'northwest','FontSize', 14,'Interpreter', 'latex')
% 
% % f=gcf;
% % exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
% %     '\Figures\about contact information vs initial condition\theta_0m10to10_theta_c-y_c_delta.png'],'Resolution',100)
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot (theta_0, y_0) vs (theta_c, y_c) in vector field (outdated version)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % let the data 'align' for reshape
% together_plot(:, together_plot(3, :) == 0.1875) = [];
% together_plot(:, together_plot(3, :) == 0.3125) = [];
% together_plot(:, together_plot(3, :) == 0.4375) = [];
% together_plot(:, together_plot(3, :) == 0.21875) = [];
% together_plot(:, together_plot(3, :) == 0.28125) = [];
% together_plot(:, together_plot(3, :) == 0.34375) = [];
% together_plot(:, together_plot(3, :) == 0.40625) = [];
% 
% % choose an 'L' to plot:
% choose_L = 0.5;
% together_plot(:, together_plot(2, :) ~= choose_L) = [];
% 
% toPlot_theta_0 = together_plot(1, :);
% toPlot_L_0 = together_plot(2, :);
% toPlot_y_0 = together_plot(3, :);
% toPlot_y_c = together_plot(5, :);
% toPlot_theta_c = together_plot(6, :);
% 
% [toPlot_y_0_X,toPlot_theta_0_Y] = meshgrid(unique(toPlot_theta_0),unique(toPlot_y_0));
% 
% toPlot_theta_c = reshape(toPlot_theta_c, [length(unique(toPlot_y_0)), length(unique(toPlot_theta_0))]);
% toPlot_y_c = reshape(toPlot_y_c, [length(unique(toPlot_y_0)), length(unique(toPlot_theta_0))]);
% toPlot_contact_X = toPlot_y_c .* sind(toPlot_theta_c);
% toPlot_contact_Y = toPlot_y_c .* cosd(toPlot_theta_c);
% 
% toPlot_contact_X_base = ~isnan(toPlot_y_c) .* ones(size(toPlot_y_c,1), size(toPlot_y_c,2)) .* sind(toPlot_theta_c);
% toPlot_contact_Y_base = ~isnan(toPlot_y_c) .* ones(size(toPlot_y_c,1), size(toPlot_y_c,2)) .* cosd(toPlot_theta_c);
% 
% figure('color', 'w'); set(gcf, 'Position', [100 100 1500 300]);
% quiver(toPlot_y_0_X,toPlot_theta_0_Y,toPlot_contact_Y_base,toPlot_contact_X_base, ...
%     'k','LineWidth', 1, 'ShowArrowHead','off'); hold on
% quiver(toPlot_y_0_X,toPlot_theta_0_Y,toPlot_contact_Y,toPlot_contact_X,'m','LineWidth', 2); axis equal
% xlabel('$Initial\ angle\ (\theta_0)$','FontSize', 18,'Interpreter', 'latex'); 
% ylabel('$Initial\ position\ (y_0)$','FontSize', 18,'Interpreter', 'latex');
% title_txt = ['$Length\ L = ',num2str(choose_L), '$'];
% title(title_txt,'FontSize', 18,'Interpreter', 'latex');
% xlim([-12 12]); ylim([-0.1 2.1]); grid on
% % f=gcf;
% % exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
% %     '\Figures\about contact information vs initial condition\contourL_', num2str(choose_L) ,'.png'],'Resolution',100)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot (theta_0, y_0) vs (theta_c, y_c) in vector field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

the_y0s = sort([0.1:0.1:0.9, 0.15:0.1:0.55, 0.325:0.05:0.525]);

toPlot_theta_0 = reshape(together_plot(1, :), 19, 8, 9);
toPlot_L_0 = reshape(together_plot(2, :), 19, 8, 9);
toPlot_y_0 = reshape(together_plot(3, :), 19, 8, 9);
toPlot_y_c = reshape(together_plot(5, :), 19, 8, 9);
toPlot_theta_c = reshape(together_plot(6, :), 19, 8, 9);

for foo = 1:numel(toPlot_y_0)
    tmp_dif = abs(toPlot_y_0(foo)-the_y0s);
    toPlot_y_0(foo) = the_y0s(tmp_dif==min(tmp_dif));
end

%%% choose one contour length to plot
% choose_L = 1; % the layer number.
% the_L_0 = unique(squeeze(toPlot_L_0(:,choose_L,:)));
% 
% toPlot_y_0_X = squeeze(toPlot_y_0(:,choose_L,:));
% toPlot_theta_0_Y = squeeze(toPlot_theta_0(:,choose_L,:));
% 
% toPlot_y_c_X = squeeze(toPlot_y_c(:,choose_L,:)) + 1; % +1: to be better shown
% toPlot_theta_c_Y = squeeze(toPlot_theta_c(:,choose_L,:));
% toPlot_contact_Y = toPlot_y_c_X .* sind(toPlot_theta_c_Y);
% toPlot_contact_X = toPlot_y_c_X .* cosd(toPlot_theta_c_Y);
% 
% 
% figure('color', 'w'); set(gcf, 'Position', [100 100 1200 800]);
% quiver(toPlot_y_0_X, toPlot_theta_0_Y/180*pi, toPlot_contact_X/50, toPlot_contact_Y/50, 0, ...
%     'm','LineWidth', 2,'ShowArrowHead','on'); 
% hold on
% axis equal;  grid on
% xlim([0.1 0.7]); ylim([-0.2 0.2]);
% ylabel('$Initial\ angle\ \theta_0\ (^\circ)$','FontSize', 18,'Interpreter', 'latex'); 
% xlabel('$Initial\ position\ y_0/h_{\rm obs}$','FontSize', 18,'Interpreter', 'latex');
% title_txt = ['$Length\ L/l_{\rm obs} = ',num2str(the_L_0), '$'];
% title(title_txt,'FontSize', 18,'Interpreter', 'latex');
% 
% ref_X = nan(size(toPlot_y_0_X, 1), size(toPlot_y_0_X, 2));
% ref_Y = ref_X;
% ref_X(16, 9) = 0.04; ref_Y(16, 9) = 0; 
% quiver(toPlot_y_0_X, toPlot_theta_0_Y/180*pi, ref_X, ref_Y, 0, 'k',...
%         'LineWidth', 2,'ShowArrowHead','off'); 
% set(gca,'TickLabelInterpreter','latex','Fontsize',14,'YTick', ...
%     -10/180*pi:2.5/180*pi:10/180*pi,'YTickLabel',-10:2.5:10);
% text(0.58,0.16, 'Reference line','Color','k','FontSize',16)
% 
% % f=gcf;
% % exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
% %     '\Figures\about contact information vs initial condition\Simu data 2023-02-23' ...
% %     '\vector map\contourL_', num2str(the_L_0) ,'.png'],'Resolution',100)


%%% plot different lengths on one figure
% legend_txt = {}; cmap = cmocean('thermal');  color_ind = 1; 
% figure('color', 'w'); set(gcf, 'Position', [100 100 1200 800]);
% for ii = 8:-1:1
% 
%     the_L_0 = unique(squeeze(toPlot_L_0(:,ii,:)));
% 
%     toPlot_y_0_X = squeeze(toPlot_y_0(:,ii,:));
%     toPlot_theta_0_Y = squeeze(toPlot_theta_0(:,ii,:));
% 
%     toPlot_y_c_X = squeeze(toPlot_y_c(:,ii,:)) + 1; % +1: to be better shown
%     toPlot_theta_c_Y = squeeze(toPlot_theta_c(:,ii,:));
%     toPlot_contact_Y = toPlot_y_c_X .* sind(toPlot_theta_c_Y);
%     toPlot_contact_X = toPlot_y_c_X .* cosd(toPlot_theta_c_Y);
% 
%     q = quiver(toPlot_y_0_X, toPlot_theta_0_Y/180*pi, toPlot_contact_X/50, toPlot_contact_Y/50, 0, 'color',...
%         cmap(color_ind*20,:),'LineWidth', 1,'ShowArrowHead','on'); 
%     q.MaxHeadSize = 0.03;
%     hold on
%     color_ind = color_ind + 1;
%     legend_txt = [legend_txt, strcat('$L/l_{\rm obs} = ',num2str(the_L_0), '$')];
% 
% end
% 
% ref_X = nan(size(toPlot_y_0_X, 1), size(toPlot_y_0_X, 2));
% ref_Y = ref_X;
% ref_X(16, 9) = 0.04; ref_Y(16, 9) = 0; 
% quiver(toPlot_y_0_X, toPlot_theta_0_Y/180*pi, ref_X, ref_Y, 0, 'k',...
%         'LineWidth', 2,'ShowArrowHead','off'); 
% hold on;
% legend_txt = [legend_txt, ''];
% axis equal;  grid on
% xlim([0.1 0.7]); ylim([-0.2 0.2]);
% ylabel('$Initial\ angle\ \theta_0\ (^\circ)$','FontSize', 18,'Interpreter', 'latex'); 
% xlabel('$Initial\ position\ y_0/h_{\rm obs}$','FontSize', 18,'Interpreter', 'latex');
% legend(legend_txt, 'Location', 'southwest','FontSize', 16,'Interpreter', 'latex')
% set(gca,'TickLabelInterpreter','latex','Fontsize',14,'YTick', ...
%     -10/180*pi:2.5/180*pi:10/180*pi,'YTickLabel',-10:2.5:10);
% text(0.58,0.16, 'Reference line','Color','k','FontSize',16)
% 
% % f=gcf;
% % exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
% %     '\Figures\about contact information vs initial condition\Simu data 2023-02-23' ...
% %     '\vector map\contourL_all.png'],'Resolution',100)


%%% plot different lengths on one figure and indicators on the arrows.
legend_txt = {}; cmap = cmocean('phase');  color_ind = 1; 
figure('color', 'w'); set(gcf, 'Position', [100 100 1200 800]);
for ii = 8:-1:1

    the_L_0 = unique(squeeze(toPlot_L_0(:,ii,:)));

    toPlot_y_0_X = squeeze(toPlot_y_0(:,ii,:));
    toPlot_theta_0_Y = squeeze(toPlot_theta_0(:,ii,:));

    toPlot_y_c_X = squeeze(toPlot_y_c(:,ii,:)) + 1; % +1: to be better shown
    toPlot_theta_c_Y = squeeze(toPlot_theta_c(:,ii,:));
    toPlot_contact_Y = toPlot_y_c_X .* sind(toPlot_theta_c_Y);
    toPlot_contact_X = toPlot_y_c_X .* cosd(toPlot_theta_c_Y);

    q = quiver(toPlot_y_0_X, toPlot_theta_0_Y/180*pi, toPlot_contact_X/50, toPlot_contact_Y/50, 0, 'color',...
        cmap(color_ind*20,:),'LineWidth', 1.5, 'ShowArrowHead','on'); 
    q.MaxHeadSize = 0.03;
    hold on
    color_ind = color_ind + 1;

    % references lines
    ref_X = toPlot_contact_X; ref_Y = toPlot_contact_Y;
    ref_X = ref_X ./ sqrt(toPlot_contact_X.^2 + toPlot_contact_Y.^2);
    ref_Y = ref_Y ./ sqrt(toPlot_contact_X.^2 + toPlot_contact_Y.^2);
    qq = quiver(toPlot_y_0_X, toPlot_theta_0_Y/180*pi, ref_X/50, ref_Y/50, 0, 'color',...
        'k','LineWidth', 0.5, 'ShowArrowHead','on');
    qq.MaxHeadSize = 0.02;

    legend_txt = [legend_txt, strcat('$L/l_{\rm obs} = ',num2str(the_L_0), '$')];
   
end

axis equal;  grid on
xlim([0.1 0.7]); ylim([-0.2 0.2]);
ylabel('$Initial\ angle\ \theta_0\ (^\circ)$','FontSize', 18,'Interpreter', 'latex'); 
xlabel('$Initial\ position\ y_0/h_{\rm obs}$','FontSize', 18,'Interpreter', 'latex');
legend_txt_tmp = legend_txt; 
legend_txt = cell(1, numel(legend_txt)*2);
for foo = 1:numel(legend_txt)
    if mod(foo, 2) == 1
        legend_txt{foo} = legend_txt_tmp{(foo+1)/2};
    else
        legend_txt{foo} = '';
    end
end
legend(legend_txt, 'Location', 'southwest','FontSize', 16,'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex','Fontsize',16,'YTick', ...
    -10/180*pi:2.5/180*pi:10/180*pi,'YTickLabel',-10:2.5:10);

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\Simu data 2023-02-23' ...
%     '\vector map\contourL_all_refontop.png'],'Resolution',100)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot (theta_0, y_0) vs (theta_c-theta_0, y_c-y_0/y_c) in vector field
% theta_c-theta_0: arrow direction
% y_c-y_0/y_c: colorcode
% Keep fiber length L = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

the_y0s = sort([0.1:0.1:0.9, 0.15:0.1:0.55, 0.325:0.05:0.525]);

toPlot_theta_0 = reshape(together_plot(1, :), 19, 8, 9);
toPlot_L_0 = reshape(together_plot(2, :), 19, 8, 9);
toPlot_y_0 = reshape(together_plot(3, :), 19, 8, 9);
toPlot_y_c = reshape(together_plot(5, :), 19, 8, 9);
% toPlot_theta_c = reshape(together_plot(6, :), 19, 8, 9);
toPlot_delta_theta = reshape(together_plot(8, :), 19, 8, 9);

for foo = 1:numel(toPlot_y_0)
    tmp_dif = abs(toPlot_y_0(foo)-the_y0s);
    toPlot_y_0(foo) = the_y0s(tmp_dif==min(tmp_dif));
end


figure('color', 'w'); set(gcf, 'Position', [100 100 1200 800]);
cmap = cmocean('balance');  

choose_L = 6; % the layer number (L = 1).
the_L_0 = unique(squeeze(toPlot_L_0(:,choose_L,:)));

toPlot_y_0_X = squeeze(toPlot_y_0(:,choose_L,:));
toPlot_theta_0_Y = squeeze(toPlot_theta_0(:,choose_L,:)); % x, y coordinates
toPlot_theta_0_Y = toPlot_theta_0_Y/180*pi;
p1 = [ toPlot_theta_0_Y(:), toPlot_y_0_X(:)]; % for arrow3 plotting

toPlot_delta_theta = squeeze(toPlot_delta_theta(:,choose_L,:));
toPlot_contact_Y = ones(size(toPlot_delta_theta, 1), ...
    size(toPlot_delta_theta, 2)) .* sind(toPlot_delta_theta);
toPlot_contact_X = ones(size(toPlot_delta_theta, 1), ...
    size(toPlot_delta_theta, 2)) .* cosd(toPlot_delta_theta); % arrow direction
p2 = [toPlot_contact_X(:), toPlot_contact_Y(:)]; % for arrow3 plotting
 
toPlot_color = squeeze(toPlot_y_c(:,choose_L,:)) - squeeze(toPlot_y_0(:,choose_L,:)); % to plot y_c-y_0
% toPlot_color = squeeze(toPlot_y_c(:,choose_L,:)); % to plot y_c 
diff_y_c = max(max(toPlot_color)) - min(min(toPlot_color));
toPlot_y_c_colorInd = round((toPlot_color - min(min(toPlot_color))) / diff_y_c * 254) + 1; 
toPlot_y_c_colorInd = toPlot_y_c_colorInd(:); % arrow color

arrow_scale = 0.06;
for ii = 1:size(p1, 1)
    if ~isnan(toPlot_y_c_colorInd(ii))
    set(gca,'ColorOrder', cmap(toPlot_y_c_colorInd(ii), :))
    arrow3(p1(ii,:),p1(ii,:)+arrow_scale*p2(ii,:), '3o', 0)
    hold on
    end
end

% hcb=colorbar('Ticks', 0:0.2:1); caxis([0 1]); cmocean('phase'); % for plot y_c
% set(hcb, 'FontSize', 24, 'TickLabelInterpreter', 'latex')
% title(hcb,'$y_{\rm c}/h_{\rm obs}$','FontSize', 24,'Interpreter', 'latex'); grid on

hcb=colorbar('Ticks', -0.4:0.2:0.4); caxis([-0.4 0.4]); cmocean('balance'); % for plot y_c-y_0
set(hcb, 'FontSize', 24, 'TickLabelInterpreter', 'latex')
title(hcb,'$\Delta{y}_{\rm c}/h_{\rm obs}$','FontSize', 24,'Interpreter', 'latex'); grid on

axis equal;  grid on; box on
ylim([0.18 0.62]); xlim([-0.25 0.3]);
xlabel('$\theta_0$','FontSize', 24,'Interpreter', 'latex'); 
ylabel('$y_0/h_{\rm obs}$','FontSize', 24,'Interpreter', 'latex');

set(gca,'TickLabelInterpreter','latex','Fontsize',24,'XTick', ...
    -10/180*pi:5/180*pi:15/180*pi,'XTickLabel',-10:5:15);

% annotation
plot([-2.5/180*pi 1/180*pi], [0.475 0.475], 'k:', 'LineWidth', 2);
annotation('textarrow',[0.575,0.469], [0.743,0.683],'String','$\Delta\theta_{\rm c}$', ...
    'FontSize', 24, 'Interpreter', 'latex')

set_plot(gcf, gca)

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\Simu data 2023-02-23' ...
%     '\vector map\vectorMap_delta-theta_yc-y0_L1_switchXY.eps'])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot (theta_0, y_0) vs theta_c in vector field and colorcode the L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

the_y0s = sort([0.1:0.1:0.9, 0.15:0.1:0.55, 0.325:0.05:0.525]);

toPlot_theta_0 = reshape(together_plot(1, :), 19, 8, 9);
toPlot_L_0 = reshape(together_plot(2, :), 19, 8, 9);
toPlot_y_0 = reshape(together_plot(3, :), 19, 8, 9);
toPlot_y_c = reshape(together_plot(5, :), 19, 8, 9);
% toPlot_theta_c = reshape(together_plot(6, :), 19, 8, 9);
toPlot_delta_theta = reshape(together_plot(8, :), 19, 8, 9);

for foo = 1:numel(toPlot_y_0)
    tmp_dif = abs(toPlot_y_0(foo)-the_y0s);
    toPlot_y_0(foo) = the_y0s(tmp_dif==min(tmp_dif));
end

%%% plot different lengths on one figure and indicators on the arrows.
legend_txt = {}; cmap = cmocean('matter');  color_ind = 1; 
figure('color', 'w'); set(gcf, 'Position', [100 100 1200 800]);
for ii = 8:-2:1

    the_L_0 = unique(squeeze(toPlot_L_0(:,ii,:)));

    toPlot_y_0_X = squeeze(toPlot_y_0(:,ii,:));
    toPlot_theta_0_Y = squeeze(toPlot_theta_0(:,ii,:));

    toPlot_theta_c_ii = squeeze(toPlot_delta_theta(:,ii,:));
    toPlot_contact_Y = ones(size(toPlot_theta_c_ii, 1), ...
        size(toPlot_theta_c_ii, 2)) .* sind(toPlot_theta_c_ii);
    toPlot_contact_X = ones(size(toPlot_theta_c_ii, 1), ...
        size(toPlot_theta_c_ii, 2)) .* cosd(toPlot_theta_c_ii); % arrow direction

    q = quiver(toPlot_theta_0_Y/180*pi, toPlot_y_0_X, toPlot_contact_X/30, toPlot_contact_Y/30, 0, 'color',...
        cmap(color_ind*60,:),'LineWidth', 2, 'ShowArrowHead','off'); 
    q.MaxHeadSize = 0.05;
    hold on
    color_ind = color_ind + 1;

    legend_txt = [legend_txt, strcat('$L/l_{\rm obs} = ',num2str(the_L_0), '$')];
   
end

% axis equal;  grid on; box on
ylim([0.08 0.62]); xlim([-0.25 0.3]);
xlabel('$\theta_0$','FontSize', 24,'Interpreter', 'latex'); 
ylabel('$y_0/h_{\rm obs}$','FontSize', 24,'Interpreter', 'latex');

% annotation
plot([0 2/180*pi], [0.475 0.475], 'k:', 'LineWidth', 2);
annotation('textarrow',[0.639,0.527], [0.773,0.727],'String','$\Delta\theta_{\rm c}$', ...
    'FontSize', 24, 'Interpreter', 'latex')

legend([legend_txt, {''}], 'Location', 'southwest','FontSize', 24,'Interpreter', 'latex')
set(gca,'TickLabelInterpreter','latex','Fontsize',24,'XTick', ...
    -10/180*pi:5/180*pi:15/180*pi,'XTickLabel',-10:5:15);

set_plot(gcf, gca)

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\Simu data 2023-02-23' ...
%     '\vector map\vectorMap_deltathetac_L_switchXY.eps'])