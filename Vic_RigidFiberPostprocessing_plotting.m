%%% Plotting for rigid fiber (after Vic_RigidFiberPostprocessing.m)
% Should be run in the path which contains the following *.mat files.

clear; close all; clc;

% basic information.
h_obs = 75; l_obs = 86.6; Obj_Mag = 0.63;

% laod data and put all in one.
load('20220624-SU8_Fibers-Individual_triangularPillar_uppoint_infos_contourL-fromAVG_with_all-Chi_timestamps.mat');
All_data1 = All_data;
load('20220913-SU8_Fibers-Individual_triangularPillar_uppoint_infos_contourL-fromAVG_with_all-Chi_timestamps.mat');
All_data2 = All_data;
load('20221004-SU8_Fibers-Individual_triangularPillar_uppoint_infos_contourL-fromAVG_with_all-Chi_timestamps.mat');
All_data3 = All_data;
load('20221005-SU8_Fibers-Individual_triangularPillar_uppoint_infos_contourL-fromAVG_with_all-Chi_timestamps.mat');
All_data4 = All_data;
All_data = [All_data1, All_data2, All_data3, All_data4];

% get names of the data.
names = [All_data(:).filename];

% get other information.
delta_y = [All_data(:).delta_y];
norm_delta_y = delta_y / h_obs; % normalized delta
contourL = [All_data(:).contour_length];
norm_contourL = contourL / l_obs; % normalized L
Chi_0 = - [All_data(:).Chi_0]; % Chi_0: minus -- to meet the definition (follow the slope of the triangular obstacle to be positive)
initial_x = [All_data(:).initial_x];
norm_initial_x = initial_x / h_obs; % normalized x_0
initial_y = [All_data(:).initial_y];
norm_initial_y = initial_y / h_obs; % normalized y_0
final_x = [All_data(:).final_x];
norm_final_x = final_x / h_obs; % normalized x_f
final_y = [All_data(:).final_y];
norm_final_y = final_y / h_obs; % normalized y_f
bypass_tip_together = [All_data(:).bypass_tip]; % if pass by the apex of the obstacle
speed = [All_data(:).ave_speed]; % average speed along X: (x_f-x_0)/(t_f-t_0)
timestamps = [All_data(:).timestamps]; % time
PoleVaulting = [All_data(:).PoleVaulting]; % if it's pole-vaulting
ApexVaulting = [All_data(:).ApexVaulting]; % if it's apex-vaulting
Sliding = [All_data(:).Sliding]; % if it's sliding
Trapping = [All_data(:).Trapping]; % if it's trapping

% calculate the delta chi (between chi_f and chi_0)
Chi = [All_data(:).Chi];
for i = 1:length(Chi)
    dt = diff(timestamps{1, i}); % [time(i+1) - time(i)]
    dchi = diff(Chi{1, i})'; % chi(i+1) - chi(i)
%     dchi_dt = dchi ./ (dt * 100); % [chi(i+1) - chi(i)] / [time(i+1) - time(i)]. Here * 100 for convenience.
    n_plus = length(find(dchi<-60));  n_minus = length(find(dchi>60)); % check if there is a rotation and determine the direction via the value of d_chi
    delta_Chi(i) = Chi{1, i}(end) - Chi{1, i}(1) + n_plus*180 - n_minus*180;
    Chi_f(i) = - Chi{1, i}(end); % the orientation of the fiber in the last frame (most downstream)
%     % plot the chi evolution.
%     figure('color', 'w'); set(gcf, 'Position', [-2000 500 2000 300]);
%     subplot(1,4,1); plot(timestamps{1, i}, Chi{1, i}); ylim([-90 90]); title('chi') % plot chi vs. time.
%     subplot(1,4,2); plot(dt); ylim([0 0.1]); title('d_t')
%     subplot(1,4,3); plot(dchi); ylim([-180 180]); title('d_c_h_i')
%     hold on; line([0 300],[60 60],'color','b');
%     hold on; line([0 300],[-60 -60],'color','b'); hold off
%     subplot(1,4,4); plot(dchi_dt); ylim([-15 15]); title('d_c_h_i/d_t')
end
delta_Chi(18) = delta_Chi(18) + 180; % correct these values.
delta_Chi(87) = delta_Chi(87)  + 180;
delta_Chi(102) = delta_Chi(102) + 180;
delta_Chi(110) = delta_Chi(110) + 180;
delta_Chi = -delta_Chi; % delta_chi: counterclockwise is positive (meet the definition).

% calculate the dx/dt (using movmean)
CoM = [All_data(:).CoM];
for j = 1:length(CoM)
    dt = diff(timestamps{1, j}); % [time(i+1) - time(i)]
    dx = diff(CoM{1, j}(1, :))'; % [x(i+1) - x(i)]
    dx_dt = movmean(dx, 7) ./ movmean(dt, 7) * Obj_Mag; % [chi(i+1) - chi(i)] / [time(i+1) - time(i)]. 
    num_up_obs = sum(CoM{1, j}(1, :) < 300); % set range (without the perturbation of the obstacle) for average speed calculation (in pixel).
    speed_upstream(j) = mean(dx_dt(1:num_up_obs)); % average speed (without the perturbation of the obstacle)
    num_down_obs = sum(CoM{1, j}(1, :) > 1749); % set range (without the perturbation of the obstacle) for average speed calculation (in pixel).
    if Trapping(j) == 1
        speed_downstream(j) = speed_upstream(j);
    else
        speed_downstream(j) = mean(dx_dt(end-num_down_obs+1:end)); % average speed (without the perturbation of the obstacle)
    end
    U0(j) = 0.5 * (speed_upstream(j) + speed_downstream(j)); % U0: average of speed_upstream and speed_downstream
    speed_ave_all(j) = mean(dx_dt); % average speed 
    interaction1(j) = 1-speed_ave_all(j)/U0(j);
    interaction2(j) = max(abs(dx_dt-U0(j)))/U0(j); % interaction2 defines as max(abs(U0-U(t)))/U0;
%     % plot the dx/dt evolution.
%     figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
%     plot(timestamps{1, j}(1:end-1), dx_dt, 'Color','m', 'LineStyle','none', 'Marker','.', 'MarkerSize', 20); 
%     xlim([0 3]); ylim([0 1200])
%     xlabel('$Time\ (s)$','FontSize', 22,'Interpreter', 'latex');
%     ylabel('$U_x\ (\mu{m}/s)$','FontSize', 22,'Interpreter', 'latex');
%     f=gcf;
%     exportgraphics(f,[num2str(j),'_',names{1, j}(1:end-4)  ,'_U_x.png'],'Resolution',100)
%     close all
end

% % % plot the dx/dt (using movmean) and chi evolution
% % % for k = 164:204
% % %     dt = diff(timestamps{1, k}); % [time(i+1) - time(i)]
% % %     dx = diff(CoM{1, k}(1, :))'; % [x(i+1) - x(i)]
% % %     dx_dt = movmean(dx, 7) ./ movmean(dt, 7) * Obj_Mag; % [chi(i+1) - chi(i)] / [time(i+1) - time(i)]. 
% % %     % plot the dx/dt evolution.
% % %     figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% % %     yyaxis left
% % %     plot(timestamps{1, k}(1:end-1), dx_dt, 'LineStyle','none', 'Marker','.', 'MarkerSize', 20); 
% % %     ylabel('$U_x\ (\mu{m}/s)$','FontSize', 22,'Interpreter', 'latex');
% % %     ylim([0 1200])
% % %     yyaxis right
% % %     plot(timestamps{1, k}(1:end-1), Chi{1, k}(1:end-1), 'LineStyle','none', 'Marker','.', 'MarkerSize', 20); 
% % %     xlim([0 3]); 
% % %     ylim([-100 100])
% % %     xlabel('$Time\ (s)$','FontSize', 22,'Interpreter', 'latex');
% % %     ylabel('$\chi\ (^{\circ})$','FontSize', 22,'Interpreter', 'latex');
% % %     f=gcf;
% % %     exportgraphics(f,['20221005_',names{1, k}(1:end-4)  ,'_U_x_Chi.png'],'Resolution',100)
% % %     close all
% % % end

% put all the information in a matrix
together = [norm_delta_y; norm_contourL; Chi_0; norm_initial_x; norm_initial_y;...
    norm_final_x; norm_final_y; bypass_tip_together; speed_ave_all;...
    delta_Chi; PoleVaulting; ApexVaulting; Sliding; Trapping; ... 
    interaction1; speed_upstream; speed_downstream; Chi_f; interaction2];
% No.9 row: ones(1, length(CoM))-speed./speed_upstream: 1 - U_bar/U0.
trapped_names = names(logical(Trapping));  % extract the trapping case names
trapped_together = together(:, logical(Trapping)); % extract the trapping case information

% select reasonable cases.
range_upstream = -5; range_downstream = 5; % set acceptable range (norm_initial_x and norm_final_x)
names(together(4, :) > range_upstream) = [];
together(:, together(4, :) > range_upstream) = []; % remove the cases too close to the obstacle on the upstream side.
names(together(6, :) < range_downstream) = [];
together(:, together(6, :) < range_downstream) = []; % remove the cases too close to the obstacle on the downstream side.

% % % %% statistics
% % % % Contour length L
% % % figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% % % edges = 20:20:180;
% % % histogram(contourL,edges);
% % % set(gca,'FontSize',16);
% % % xlabel('$Contour\ length\ L\ ({\mu}m)$','FontSize', 22,'Interpreter', 'latex');
% % % ylabel('$Number\ of\ cases$','FontSize', 22,'Interpreter', 'latex');
% % % xlim([20 160]);
% % % % f=gcf;
% % % % exportgraphics(f,'Statistics_contourL.png','Resolution',100)
% % % 
% % % % Initial angle chi_0
% % % figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% % % edges = -90:20:90;
% % % histogram(Chi_0,edges);
% % % set(gca,'FontSize',16);
% % % xlabel('$Initial\ angle\ \chi_0\ (^{\circ})$','FontSize', 22,'Interpreter', 'latex');
% % % ylabel('$Number\ of\ cases$','FontSize', 22,'Interpreter', 'latex');
% % % xlim([-90 90]);
% % % % f=gcf;
% % % % exportgraphics(f,'Statistics_Chi0.png','Resolution',100)
% % % 
% % % % Normalized initial position (y_0/h_obs)
% % % figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% % % edges = -0.2:0.2:1.2;
% % % histogram(norm_initial_y,edges);
% % % set(gca,'FontSize',16);
% % % xlabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% % % ylabel('$Number\ of\ cases$','FontSize', 22,'Interpreter', 'latex');
% % % xlim([-0.2 1.2]);
% % % % f=gcf;
% % % % exportgraphics(f,'Statistics_y0.png','Resolution',100)



%% plotting
% create dialog box to gather user input for plotting;
names_plot = names; together_plot = together;
prompt = {'The lower bound of the initial angle:', 'The upper bound of the initial angle:', ...
    'The lower bound of the contour length:','The upper bound of the contour length:'...
    'The lower bound of the initial position:','The upper bound of the initial position:'};
definput = {'-10', '10', '0.5', '1', '0', '1'};
answer = inputdlg(prompt, 'Input (please input NaN if there is no bound)', [1 35] , definput);

%%%%%%%%%%%%%%%%%%%%%%%%% assign the values %%%%%%%%%%%%%%%%%%%%%%%%%%
range_chi0_low = str2double(answer{1,1}); if isnan(range_chi0_low); range_chi0_low = -91; end
range_chi0_up = str2double(answer{2,1}); if isnan(range_chi0_up); range_chi0_up = 91; end
range_L_low = str2double(answer{3,1}); if isnan(range_L_low); range_L_low = 0; end
range_L_up = str2double(answer{4,1}); if isnan(range_L_up); range_L_up = 10; end
range_y0_low = str2double(answer{5,1}); if isnan(range_y0_low); range_y0_low = -10; end
range_y0_up = str2double(answer{6,1}); if isnan(range_y0_up); range_y0_up = -91; end
% chi_0
names_plot(together_plot(3, :) < range_chi0_low) = [];
together_plot(:, together_plot(3, :) < range_chi0_low) = []; 
names_plot(together_plot(3, :) > range_chi0_up) = [];
together_plot(:, together_plot(3, :) > range_chi0_up) = [];
% L
names_plot(together_plot(2, :) < range_L_low) = [];
together_plot(:, together_plot(2, :) < range_L_low) = []; % the contour length: 0.5 < L < 1
names_plot(together_plot(2, :) > range_L_up) = [];
together_plot(:, together_plot(2, :) > range_L_up) = [];
% y_0
names_plot(together_plot(5, :) < 0) = [];
together_plot(:, together_plot(5, :) < 0) = [];  % the initial position: 0 < y_0 < 1
names_plot(together_plot(5, :) > 1) = [];
together_plot(:, together_plot(5, :) > 1) = [];

%%%%%%%%%%%%%%%%%%%%%%%%% scatter plots %%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the deviation vs L & y_0 (without classification):
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
cmap = cmocean('thermal');
scatter(together_plot(2, :), together_plot(5, :), 150, together_plot(1, :), 'Filled')
cmap(size(together_plot,2)); 
hcb=colorbar;
title(hcb,'$Deviation\ (\delta/h_{obs})$','FontSize', 16,'Interpreter', 'latex'); grid on
set(gca,'FontSize',16);
xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
if range_L_low ~= 0; xlim([range_L_low range_L_up]); end
if range_y0_low ~= -10; ylim([range_y0_low range_y0_up]); end
% f=gcf;
% exportgraphics(f,'delta_vs_L-y0.png','Resolution',100)
% % %%%%%% pick the data point on the last plot to get the case name %%%%%%%
% % [the_L, the_y0] = ginput(1); % pick up the point you want to show the trajectory.
% % the_loc = intersect(find(together_plot(2, :)>the_L*0.98 & together_plot(2, :)<the_L*1.02),...
% %     find(together_plot(5, :)>the_y0*0.98 & together_plot(5, :)<the_y0*1.02));  % Change the control range if there is an error.
% % names_plot(the_loc)     % show the file

% plot the deviation vs L & y_0 (with classification):
figure('color', 'w'); set(gcf, 'Position', [100 100 1200 600]);
cmap = cmocean('thermal');
bypass_edge_together = together_plot(:, ~logical(together_plot(8, :)));
bypass_tip_together = together_plot(:, logical(together_plot(8, :)));
pole_vaulting_together = together_plot(:, logical(together_plot(11, :))); 
apex_vaulting_together = together_plot(:, logical(together_plot(12, :)));
sliding_together = together_plot(:, logical(together_plot(13, :)));
scatter(trapped_together(2, :), trapped_together(5, :), 150, trapped_together(1, :), 'Filled', 'diamond'); hold on 
scatter(bypass_edge_together(2, :), bypass_edge_together(5, :), 150, bypass_edge_together(1, :), 'Filled'); hold on
scatter(bypass_tip_together(2, :), bypass_tip_together(5, :), 150, bypass_tip_together(1, :), 'Filled', 'square'); hold on
scatter(pole_vaulting_together(2, :), pole_vaulting_together(5, :), 150, pole_vaulting_together(1, :), 'Filled', '^'); hold on
scatter(apex_vaulting_together(2, :), apex_vaulting_together(5, :), 150, apex_vaulting_together(1, :), 'Filled', 'v'); hold on
scatter(sliding_together(2, :), sliding_together(5, :), 150, sliding_together(1, :), 'Filled', 'pentagram');
cmap(size(together_plot,2)); 
hcb=colorbar;
title(hcb,'$Deviation\ (\delta/h_{obs})$','FontSize', 16,'Interpreter', 'latex'); grid on
set(gca,'FontSize',16);
xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
if range_L_low ~= 0; xlim([range_L_low range_L_up]); end
if range_y0_low ~= -10; ylim([range_y0_low range_y0_up]); end
% f=gcf;
% exportgraphics(f,'delta_vs_L-y0_classification.png','Resolution',100)

% plot the y_0 vs deviation & (1-U_bar/U0) <interaction1> (with classification):
figure('color', 'w'); set(gcf, 'Position', [100 100 1200 600]);
cmap = cmocean('thermal');
bypass_edge_together = together_plot(:, ~logical(together_plot(8, :)));
bypass_tip_together = together_plot(:, logical(together_plot(8, :)));
pole_vaulting_together = together_plot(:, logical(together_plot(11, :))); 
apex_vaulting_together = together_plot(:, logical(together_plot(12, :)));
sliding_together = together_plot(:, logical(together_plot(13, :)));
scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond'); hold on  % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'square'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', '^'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'v'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'pentagram'); hold on % for legend only
scatter(trapped_together(15, :), trapped_together(1, :), 200, trapped_together(5, :), 'Filled', 'diamond'); hold on 
scatter(bypass_edge_together(15, :), bypass_edge_together(1, :), 200, bypass_edge_together(5, :), 'Filled'); hold on
scatter(bypass_tip_together(15, :), bypass_tip_together(1, :), 200, bypass_tip_together(5, :), 'Filled', 'square'); hold on
scatter(pole_vaulting_together(15, :), pole_vaulting_together(1, :), 200, pole_vaulting_together(5, :), 'Filled', '^'); hold on
scatter(apex_vaulting_together(15, :), apex_vaulting_together(1, :), 200, apex_vaulting_together(5, :), 'Filled', 'v'); hold on
scatter(sliding_together(15, :), sliding_together(1, :), 200, sliding_together(5, :), 'Filled', 'pentagram');
cmap(size(together_plot,2)); 
hcb=colorbar;
title(hcb,'$Initial\ position\ (y_0/h_{obs})$','FontSize', 16,'Interpreter', 'latex'); grid on
set(gca,'FontSize',16);
xlabel('$1-\bar{U}/U_0$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Deviation\ (\delta/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
legend({'Trapping','Passing Below','Passing Above (apex side)','Pole Vaulting','Apex Vaulting', ...
    'Sliding'}, 'Location', 'southwest','FontSize', 14,'Interpreter', 'latex')
% if range_L_low ~= 0; xlim([range_L_low range_L_up]); end
% if range_y0_low ~= -10; ylim([range_y0_low range_y0_up]); end
% f=gcf;
% exportgraphics(f,'y0_vs_interaction-delta_classification.png','Resolution',100)

% plot the y_0 vs deviation & (1-U_bar/U0) <interaction1> (with classification and further data cleaning !!!):
figure('color', 'w'); set(gcf, 'Position', [100 100 1200 600]);
cmap = cmocean('thermal');
together_plot_filtered = together_plot;
together_plot_filtered(:, together_plot_filtered(9, :) > 800) = []; % remove the cases that are too fast
together_plot_filtered(:, together_plot_filtered(9, :) < 400) = []; % remove the cases that are too slow
abs_delta_U = abs(together_plot_filtered(17, :) - together_plot_filtered(16, :)); % |U_f - u_0|
abs_delta_chi = abs(together_plot_filtered(18, :) - together_plot_filtered(3, :)); % |chi_f - chi_0|
together_plot_filtered(:, and(abs_delta_U > 100, abs_delta_chi < 10)) = []; % remove the cases that have large U-differences but very small chi-differences.
bypass_edge_together = together_plot_filtered(:, ~logical(together_plot_filtered(8, :)));
bypass_tip_together = together_plot_filtered(:, logical(together_plot_filtered(8, :)));
pole_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(11, :))); 
apex_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(12, :)));
sliding_together = together_plot_filtered(:, logical(together_plot_filtered(13, :)));
scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond'); hold on  % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'square'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', '^'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'v'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'pentagram'); hold on % for legend only
scatter(trapped_together(15, :), trapped_together(1, :), 200, trapped_together(5, :), 'Filled', 'diamond'); hold on 
scatter(bypass_edge_together(15, :), bypass_edge_together(1, :), 200, bypass_edge_together(5, :), 'Filled'); hold on
scatter(bypass_tip_together(15, :), bypass_tip_together(1, :), 200, bypass_tip_together(5, :), 'Filled', 'square'); hold on
scatter(pole_vaulting_together(15, :), pole_vaulting_together(1, :), 200, pole_vaulting_together(5, :), 'Filled', '^'); hold on
scatter(apex_vaulting_together(15, :), apex_vaulting_together(1, :), 200, apex_vaulting_together(5, :), 'Filled', 'v'); hold on
scatter(sliding_together(15, :), sliding_together(1, :), 200, sliding_together(5, :), 'Filled', 'pentagram');
cmap(size(together_plot_filtered,2)); 
hcb=colorbar; caxis([0 1])
title(hcb,'$Initial\ position\ (y_0/h_{obs})$','FontSize', 16,'Interpreter', 'latex'); grid on
set(gca,'FontSize',16);
xlabel('$1-\bar{U}/U_0$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Deviation\ (\delta/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
legend({'Trapping','Passing Below','Passing Above (apex side)','Pole Vaulting','Apex Vaulting', ...
    'Sliding'}, 'Location', 'southwest','FontSize', 14,'Interpreter', 'latex')
% if range_L_low ~= 0; xlim([range_L_low range_L_up]); end
% if range_y0_low ~= -10; ylim([range_y0_low range_y0_up]); end
% f=gcf;
% exportgraphics(f,'y0_vs_interaction-delta_classification_datacleaning.png','Resolution',100)

% plot the y_0 vs deviation & max(abs(U0-U(t)))/U0 <interaction2> (with classification and further data cleaning !!!):
figure('color', 'w'); set(gcf, 'Position', [100 100 1200 600]);
% cmap = cmocean('thermal');
cmap = colormap("jet");
together_plot_filtered = together_plot;
together_plot_filtered(:, together_plot_filtered(9, :) > 800) = []; % remove the cases that are too fast
together_plot_filtered(:, together_plot_filtered(9, :) < 400) = []; % remove the cases that are too slow
abs_delta_U = abs(together_plot_filtered(17, :) - together_plot_filtered(16, :)); % |U_f - u_0|
abs_delta_chi = abs(together_plot_filtered(18, :) - together_plot_filtered(3, :)); % |chi_f - chi_0|
together_plot_filtered(:, and(abs_delta_U > 100, abs_delta_chi < 10)) = []; % remove the cases that have large U-differences but very small chi-differences.
bypass_edge_together = together_plot_filtered(:, ~logical(together_plot_filtered(8, :))); 
bypass_tip_together = together_plot_filtered(:, logical(together_plot_filtered(8, :)));
pole_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(11, :))); 
apex_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(12, :)));
sliding_together = together_plot_filtered(:, logical(together_plot_filtered(13, :)));
scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond'); hold on  % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'square'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', '^'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', 'v'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', 'pentagram'); hold on % for legend only
scatter(trapped_together(19, :), trapped_together(1, :), 200, trapped_together(5, :), 'Filled', 'diamond'); hold on 
scatter(bypass_edge_together(19, :), bypass_edge_together(1, :), 200, bypass_edge_together(5, :), 'Filled'); hold on
scatter(bypass_tip_together(19, :), bypass_tip_together(1, :), 200, bypass_tip_together(5, :), 'Filled', 'square'); hold on
scatter(pole_vaulting_together(19, :), pole_vaulting_together(1, :), 200, pole_vaulting_together(5, :), 'Filled', '^'); hold on
% scatter(apex_vaulting_together(19, :), apex_vaulting_together(1, :), 200, apex_vaulting_together(5, :), 'Filled', 'v'); hold on
% scatter(sliding_together(19, :), sliding_together(1, :), 200, sliding_together(5, :), 'Filled', 'pentagram');
cmap(size(together_plot_filtered,2)); 
hcb=colorbar; caxis([0 1])
title(hcb,'$Initial\ position\ (y_0/h_{obs})$','FontSize', 16,'Interpreter', 'latex'); grid on
set(gca,'FontSize',16);
xlabel('$max\left|U_0-U(t)\right|/U_0$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Deviation\ (\delta/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
legend({'Trapping','Below','Above','Pole-vaulting'...
    }, 'Location', 'southwest','FontSize', 14,'Interpreter', 'latex')
% title('$-10<{\chi}_0<10,\ 0.5<L<1$','FontSize', 22,'Interpreter', 'latex')
% if range_L_low ~= 0; xlim([range_L_low range_L_up]); end
% if range_y0_low ~= -10; ylim([range_y0_low range_y0_up]); end
% f=gcf;
% exportgraphics(f,'y0_vs_interaction2-delta_classification_datacleaning_Chi0-m10to10_L-0p5to1.png','Resolution',100)
[the_inter2, the_delta] = ginput(1); % pick up the point you want to show the trajectory.
the_loc = intersect(find(together_plot(19, :)>the_inter2*0.98 & together_plot(19, :)<the_inter2*1.02),...
    find(together_plot(1, :)>the_delta*0.98 & together_plot(1, :)<the_delta*1.02));  % Change the control range if there is an error.
names_plot(the_loc)     % show the file

% plot the speed vs L & y_0: 
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
cmap = cmocean('thermal');
scatter(together_plot(2, :), together_plot(5, :), 200, together_plot(9, :), 'Filled', 'hexagram') % remember to change the 9th row of 'together'
cmap(size(together_plot,2)); 
hcb=colorbar;
title(hcb,'$Average\ speed\ (\mu{m}/s)$','FontSize', 16,'Interpreter', 'latex');
grid on
set(gca,'FontSize',16);
xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,'Ubar_vs_L-y0.png','Resolution',100)

% plot (1-U_bar/U0) <interaction1> vs L & y_0:
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
cmap = cmocean('thermal');
scatter(together_plot(2, :), together_plot(5, :), 200, together_plot(15, :), 'Filled', 'hexagram') % remember to change the 9th row of 'together'
cmap(size(together_plot,2)); 
hcb=colorbar;
title(hcb,'$1-\bar{U}/U_0$','FontSize', 16,'Interpreter', 'latex');
grid on
set(gca,'FontSize',16);
xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,'interaction_vs_L-y0.png','Resolution',100)

% plot the delta Chi vs L & y_0:
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
cmap = cmocean('thermal');
scatter(together_plot(2, :), together_plot(5, :), 200, together_plot(10, :)', 'Filled', 'hexagram')    % add ' to avoid obfuscating to a triplet.
cmap(size(together_plot,2)); 
hcb=colorbar;
title(hcb,'$\chi_f - \chi_0\ (^{\circ})$','FontSize', 16,'Interpreter', 'latex');
grid on
set(gca,'FontSize',16);
xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,'deltaChi_vs_L-y0.png','Resolution',100)



%%%%%%%%%%%%%%%%%%%%%%%%% other plots %%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the |delta Chi| vs deviation 
figure('color', 'w'); set(gcf, 'Position', [100 100 600 400]);
plot(abs(together_plot(10, :)), together_plot(1, :), 'Color','r', 'LineStyle','none', 'Marker','*', 'MarkerSize', 5)
set(gca,'FontSize',16);
xlabel('$\left|\chi_f - \chi_0 \right| (^{\circ})$','FontSize',22,'Interpreter', 'latex');
ylabel('$Deviation\ (\delta/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,'absDeltaChi_vs_delta.png','Resolution',100)

% plot the (1 - U_bar/U0) <interaction1>  vs deviation
figure('color', 'w'); set(gcf, 'Position', [100 100 600 400]);
plot(together_plot(15, :), together_plot(1, :), 'Color','r', 'LineStyle','none', 'Marker','*', 'MarkerSize', 5)
set(gca,'FontSize',16);
xlabel('$1-\bar{U}/U_0$','FontSize',22,'Interpreter', 'latex');
ylabel('$Deviation\ (\delta/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,'interaction_vs_delta.png','Resolution',100)

% plot the (speed_downstream-speed_upstream) vs Chi_0/Chi_f
figure('color', 'w'); set(gcf, 'Position', [100 100 600 400]);
yyaxis left
plot((together_plot(17, :) - together_plot(16, :)), together_plot(3, :), 'LineStyle','none', 'Marker','*', 'MarkerSize', 7)
ylabel('$\chi_0\ (^{\circ})$','FontSize', 22,'Interpreter', 'latex');
yyaxis right
plot((together_plot(17, :) - together_plot(16, :)), together_plot(18, :), 'LineStyle','none', 'Marker','o', 'MarkerSize', 7)
ylabel('$\chi_f\ (^{\circ})$','FontSize', 22,'Interpreter', 'latex');
set(gca,'FontSize',16);
xlabel('$U_{downstream}-U_{upstream}\ (\mu{m}/s)$','FontSize',22,'Interpreter', 'latex');
xlim([-300 300])
% f=gcf;
% exportgraphics(f,'deltaU_vs_chi0__eltaU_vs_chif.png','Resolution',100)