%%% Plotting for rigid fiber (after Vic_RigidFiberPostprocessing.m)
% Should be run in the path which contains the following *.mat files.

clear; close all; clc;

% basic information.
h_obs = 75; l_obs = 86.6; Obj_Mag = 0.63;

% laod data and put all in one.
load('F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\20220624-SU8_Fibers-Individual_triangularPillar_uppoint_information_full.mat');
All_data1 = delta_correction(All_data, 0);
load('F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\20220913-SU8_Fibers-Individual_triangularPillar_uppoint_information_full.mat');
All_data2 = All_data;
load('F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\20221004-SU8_Fibers-Individual_triangularPillar_uppoint_information_full.mat');
All_data3 = All_data;
load('F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\20221005-SU8_Fibers-Individual_triangularPillar_uppoint_information_full.mat');
All_data4 = All_data;
load('F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\20230102-SU8_Fibers-Individual_triangularPillar_uppoint_information_full.mat');
All_data5 = delta_correction(All_data, 0);
load('F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\202208-202209_SU8fibre-and-Tracer-Individual_triangularPillar_uppoint_information_full.mat');
All_data6 = All_data;
All_data = [All_data1, All_data2, All_data3, All_data4, All_data5, All_data6];

% get names of the data.
names = [All_data(:).filename];
acceptability = cellfun(@isnumeric, [All_data(:).acceptability]);

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
if_bypass_tip_together = [All_data(:).bypass_tip]; % if pass by the apex of the obstacle
ave_speed = [All_data(:).ave_speed]; % average speed along X: (x_f-x_0)/(t_f-t_0)
speed_ave_all = [All_data(:).speed_ave_all]; % average speed along X: mean(dx/dt) 
speed_downstream = [All_data(:).speed_downstream];
speed_upstream = [All_data(:).speed_upstream];
timestamps = [All_data(:).timestamps]; % time
PoleVaulting = [All_data(:).PoleVaulting]; % if it's pole-vaulting
ApexVaulting = [All_data(:).ApexVaulting]; % if it's apex-vaulting
Sliding = [All_data(:).Sliding]; % if it's sliding
Trapping = [All_data(:).Trapping]; % if it's trapping
interaction1 = [All_data(:).interaction1]; % interaction index 1
interaction2 = [All_data(:).interaction2]; % interaction index 2
interaction3 = [All_data(:).interaction3]; % interaction index 3
CoMs = [All_data(:).CoM];
Foo = [All_data(:).y_c]+10; Foo(isnan(Foo)) = 0;
if_contact = logical(Foo); % if fiber contacts the obstacle

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
% % % CoM = [All_data(:).CoM];
% % % for j = 1:length(CoM)
% % %     dt = diff(timestamps{1, j}); % [time(i+1) - time(i)]
% % %     dx = diff(CoM{1, j}(1, :))'; % [x(i+1) - x(i)]
% % %     dx_dt = movmean(dx, 7) ./ movmean(dt, 7) * Obj_Mag; % [chi(i+1) - chi(i)] / [time(i+1) - time(i)]. 
% % % %     num_up_obs = sum(CoM{1, j}(1, :) < 300); % set range (without the perturbation of the obstacle) for average speed calculation (in pixel).
% % % %     speed_upstream(j) = mean(dx_dt(1:num_up_obs)); % average speed (without the perturbation of the obstacle)
% % % %     num_down_obs = sum(CoM{1, j}(1, :) > 1749); % set range (without the perturbation of the obstacle) for average speed calculation (in pixel).
% % % %     if Trapping(j) == 1
% % % %         speed_downstream(j) = speed_upstream(j);
% % % %     else
% % % %         speed_downstream(j) = mean(dx_dt(end-num_down_obs+1:end)); % average speed (without the perturbation of the obstacle)
% % % %     end
% % % %     U0(j) = 0.5 * (speed_upstream(j) + speed_downstream(j)); % U0: average of speed_upstream and speed_downstream
% % % %     speed_ave_all(j) = mean(dx_dt); % average speed 
% % % %     interaction1(j) = 1-speed_ave_all(j)/U0(j);
% % % %     interaction2(j) = max(abs(dx_dt-U0(j)))/U0(j); % interaction2 defines as max(abs(U0-U(t)))/U0;
% % %     % plot the dx/dt evolution.
% % %     figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% % %     plot(timestamps{1, j}(1:end-1), dx_dt, 'Color','m', 'LineStyle','none', 'Marker','.', 'MarkerSize', 20); 
% % %     xlim([0 3]); ylim([0 1200])
% % %     xlabel('$Time\ (s)$','FontSize', 22,'Interpreter', 'latex');
% % %     ylabel('$U_x\ (\mu{m}/s)$','FontSize', 22,'Interpreter', 'latex');
% % % %     f=gcf;
% % % %     exportgraphics(f,[num2str(j),'_',names{1, j}(1:end-4)  ,'_U_x.png'],'Resolution',100)
% % %     close all
% % % end

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
% % %     ylabel('$\theta\ (^{\circ})$','FontSize', 22,'Interpreter', 'latex');
% % %     f=gcf;
% % %     exportgraphics(f,['20221005_',names{1, k}(1:end-4)  ,'_U_x_Chi.png'],'Resolution',100)
% % %     close all
% % % end

% put all the information in a matrix
together = [norm_delta_y; norm_contourL; Chi_0; norm_initial_x; norm_initial_y;...
    norm_final_x; norm_final_y; if_bypass_tip_together; ave_speed; speed_ave_all; ...
    delta_Chi; PoleVaulting; ApexVaulting; Sliding; Trapping; ... 
    Chi_f; interaction1; interaction2; interaction3; speed_upstream; speed_downstream; if_contact];
% %     No.1 row: normalized delta
% %     No.2 row: normalized L
% %     No.3 row: chi_0
% %     No.4 row: normalized x_0
% %     No.5 row: normalized y_0
% %     No.6 row: normalized x_f
% %     No.7 row: normalized y_f
% %     No.8 row: if pass by the apex of the obstacle
% %     No.9 row: average speed along X: (x_f-x_0)/(t_f-t_0) 
% %     No.10 row: average speed along X: mean(dx/dt) 
% %     No.11 row: delta chi
% %     No.12 row: if it's pole-vaulting
% %     No.13 row: if it's apex-vaulting
% %     No.14 row: if it's sliding
% %     No.15 row: if it's trapping
% %     No.16 row: chi_f
% %     No.17 row: interaction index 1
% %     No.18 row: interaction index 2
% %     No.19 row: interaction index 3
% %     No.20 row: upstream speed
% %     No.21 row: downstream speed
% %     No.22 row: if fiber contacts the obstacle

names(~logical(acceptability)) = [];
together(:, ~logical(acceptability)) = [];
Trapping(~logical(acceptability)) = [];
% these cases are not in the channel mid-plane by manual selection.
% In the excel of PoleVaultingCheck.xls, their 'Deviation correction based on tracers' are NaN.

trapped_names = names(logical(Trapping));  % extract the trapping case names
trapped_together = together(:, logical(Trapping)); % extract the trapping case information
names(logical(Trapping)) = [];
together(:, logical(Trapping)) = []; % remove trapping cases

% select reasonable cases (use this loop only for interaction index 3).
% ifOK = zeros(1, length(CoMs));
% for i = 1:length(CoMs)
%     xxx = CoMs{1, i}(1, :);
%     if i < 26
%         num = length(find(xxx > 1030 & xxx < 1170));  % the numbers are the CoM of the obstacle.
%     elseif i < 111
%         num = length(find(xxx > 993 & xxx < 1133));
%     elseif i < 164
%         num = length(find(xxx > 922 & xxx < 1062));
%     else
%         num = length(find(xxx > 934 & xxx < 1074));
%     end
% 
%     if num > 12  % the acceptable minimum frame number in the near-obstacle-region.
%         ifOK(i) = 1;
%     end
% end

% names = names(logical(ifOK));  
% together = together(:, logical(ifOK));

% select reasonable cases.
range_upstream = -3; range_downstream = 3; % set acceptable range (norm_initial_x and norm_final_x)
names(together(4, :) > range_upstream) = [];
together(:, together(4, :) > range_upstream) = []; % remove the cases too close to the obstacle on the upstream side.
names(together(6, :) < range_downstream) = [];
together(:, together(6, :) < range_downstream) = []; % remove the cases too close to the obstacle on the downstream side.

% % % %% statistics
% % % % Contour length L
% % % figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% % % edges = [0:17.32:190.52]/l_obs;
% % % histogram(contourL/l_obs,edges);
% % % set(gca,'FontSize',16);
% % % xlabel('$Contour\ length\ L$','FontSize', 22,'Interpreter', 'latex');
% % % ylabel('$Number\ of\ cases$','FontSize', 22,'Interpreter', 'latex');
% % % xlim([0 2.2]);
% % % % f=gcf;
% % % % exportgraphics(f,'Statistics_contourL.png','Resolution',100)
% % % 
% % % % Initial angle chi_0
% % % figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% % % edges = -90:20:90;
% % % histogram(Chi_0,edges);
% % % set(gca,'FontSize',16);
% % % xlabel('$Initial\ angle\ \theta_0\ (^{\circ})$','FontSize', 22,'Interpreter', 'latex');
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

% create dialog box to gather user input for plotting;
names_plot = names; together_plot = together;
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
range_L_up = str2double(answer{4,1}); if isnan(range_L_up); range_L_up = 10; end
range_y0_low = str2double(answer{5,1}); if isnan(range_y0_low); range_y0_low = -10; end
range_y0_up = str2double(answer{6,1}); if isnan(range_y0_up); range_y0_up = 10; end
% chi_0
names_plot(together_plot(3, :) < range_chi0_low) = [];
together_plot(:, together_plot(3, :) < range_chi0_low) = []; 
names_plot(together_plot(3, :) > range_chi0_up) = [];
together_plot(:, together_plot(3, :) > range_chi0_up) = [];
% L
names_plot(together_plot(2, :) < range_L_low) = [];
together_plot(:, together_plot(2, :) < range_L_low) = []; 
names_plot(together_plot(2, :) > range_L_up) = [];
together_plot(:, together_plot(2, :) > range_L_up) = [];
% y_0
names_plot(together_plot(5, :) < range_y0_low) = [];
together_plot(:, together_plot(5, :) < range_y0_low) = [];  
names_plot(together_plot(5, :) > range_y0_up) = [];
together_plot(:, together_plot(5, :) > range_y0_up) = [];
% chi_0 for trapping cases
trapped_names(trapped_together(3, :) < range_chi0_low) = [];
trapped_together(:, trapped_together(3, :) < range_chi0_low) = []; 
trapped_names(trapped_together(3, :) > range_chi0_up) = [];
trapped_together(:, trapped_together(3, :) > range_chi0_up) = [];
% L for trapping cases
trapped_names(trapped_together(2, :) < range_L_low) = [];
trapped_together(:, trapped_together(2, :) < range_L_low) = []; 
trapped_names(trapped_together(2, :) > range_L_up) = [];
trapped_together(:, trapped_together(2, :) > range_L_up) = [];
% y_0 for trapping cases
trapped_names(trapped_together(5, :) < range_y0_low) = [];
trapped_together(:, trapped_together(5, :) < range_y0_low) = [];  
trapped_names(trapped_together(5, :) > range_y0_up) = [];
trapped_together(:, trapped_together(5, :) > range_y0_up) = [];

%% Save data to AllinONE (with further data cleaning) !!!!!!!!!!!!!!!!!!!!!!!:  

together_plot_filtered = together_plot;
names_plot_filtered = names_plot;
% together_plot_filtered(:, together_plot_filtered(10, :) > 800) = []; % remove the cases that are too fast
% together_plot_filtered(:, together_plot_filtered(10, :) < 400) = []; % remove the cases that are too slow
abs_delta_U = abs(together_plot_filtered(21, :) - together_plot_filtered(20, :)); % |U_f - u_0|
abs_delta_chi = abs(together_plot_filtered(16, :) - together_plot_filtered(3, :)); % |chi_f - chi_0|

together_plot_filtered(:, and(abs_delta_U > 100, abs_delta_chi < 10)) = []; % remove the cases that have large U-differences but very small chi-differences.
names_plot_filtered(:, and(abs_delta_U > 100, abs_delta_chi < 10)) = [];

% % %% save the data
% % Trappping = [trapped_together(3, :)', trapped_together(5, :)', trapped_together(1, :)', trapped_together(2, :)'];
% % Below = [bypass_edge_together(3, :)', bypass_edge_together(5, :)', bypass_edge_together(1, :)', bypass_edge_together(2, :)'];
% % Above = [bypass_tip_together(3, :)', bypass_tip_together(5, :)', bypass_tip_together(1, :)', bypass_tip_together(2, :)'];
% % Pole_vaulting = [pole_vaulting_together(3, :)', pole_vaulting_together(5, :)', pole_vaulting_together(1, :)', pole_vaulting_together(2, :)'];
% % readme = '1st column: theta_0;   2nd col: y_0;   3rd col: delta;   4th col: L.';
% % save('D:\Dropbox\Collaboration - LadHyX\Give_to_Clement\Data & Figures\delta_vs_theta0-y0_alldata_20230110.mat','Trappping','Below','Above','Pole_vaulting', 'readme');

% % %% save the data (All in one)
AllinONE.name = [names_plot_filtered, trapped_names];
AllinONE.delta = [together_plot_filtered(1, :), trapped_together(1, :)]; % normalized delta by h_obs
AllinONE.L = [together_plot_filtered(2, :), trapped_together(2, :)]; % normalized L by l_obs
AllinONE.theta_0 = [together_plot_filtered(3, :), trapped_together(3, :)]; % theta_0
AllinONE.x_0 = [together_plot_filtered(4, :), trapped_together(4, :)]; % normalized x_0 by h_obs (initial position)
AllinONE.y_0 = [together_plot_filtered(5, :), trapped_together(5, :)]; % normalized y_0 by h_obs
AllinONE.x_f = [together_plot_filtered(6, :), trapped_together(6, :)]; % normalized x_f by h_obs (final position)
AllinONE.y_f = [together_plot_filtered(7, :), trapped_together(7, :)]; % normalized y_f by h_obs
AllinONE.if_above = [together_plot_filtered(8, :), trapped_together(8, :)]; % if 'Above'
AllinONE.if_PoleVaulting = [together_plot_filtered(12, :), trapped_together(12, :)]; % if 'PoleVaulting'
AllinONE.if_Trapping = [together_plot_filtered(15, :), trapped_together(15, :)]; % if 'Trapping'
AllinONE.if_Contacting = [together_plot_filtered(22, :), trapped_together(22, :)]; % if 'Contacting'

readme = 'See struct AllinONE, field name gives the data meaning!';
save('D:\Dropbox\Collaboration - LadHyX\Give_to_Clement\FSI - Rigid Fiber &  Individual Obstacle\Rigid_Fiber_expdata_20230629.mat','AllinONE', 'readme');


%% plotting
%%%%%%%%%%%%%%%%%%%%%%%%% scatter plots %%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the deviation vs L & y_0 (without classification):
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% cmap = cmocean('thermal');
cmap = colormap("jet");
scatter(together_plot(2, :)', together_plot(5, :)', 100, together_plot(1, :)', 'Filled', 'MarkerEdgeColor','k')
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
%%%%%% pick the data point on the last plot to get the case name %%%%%%%
% % [the_L, the_y0] = ginput(1); % pick up the point you want to show the trajectory.
% % the_loc = intersect(find(together_plot(2, :)>the_L*0.98 & together_plot(2, :)<the_L*1.02),...
% %     find(together_plot(5, :)>the_y0*0.98 & together_plot(5, :)<the_y0*1.02));  % Change the control range if there is an error.
% % names_plot(the_loc)     % show the file

%% plot the deviation vs L & y_0 (with classification and further data cleaning (optional)):
figure('color', 'w'); set(gcf, 'Position', [100 100 1000 500]);
% cmap = cmocean('thermal');
cmap = colormap("jet");
together_plot_filtered = together_plot;
% together_plot_filtered(:, together_plot_filtered(10, :) > 800) = []; % remove the cases that are too fast
% together_plot_filtered(:, together_plot_filtered(10, :) < 400) = []; % remove the cases that are too slow
abs_delta_U = abs(together_plot_filtered(21, :) - together_plot_filtered(20, :)); % |U_f - u_0|
abs_delta_chi = abs(together_plot_filtered(16, :) - together_plot_filtered(3, :)); % |chi_f - chi_0|
together_plot_filtered(:, and(abs_delta_U > 100, abs_delta_chi < 10)) = []; % remove the cases that have large U-differences but very small chi-differences.
onlypass_ind = ~(logical(together_plot_filtered(12, :))); 
% the case index without pole-vaulting, apex-vaulting and sliding.
bypass_edge_together = together_plot_filtered(:, and(~logical(together_plot_filtered(8, :)), onlypass_ind));
bypass_tip_together = together_plot_filtered(:, and(logical(together_plot_filtered(8, :)), onlypass_ind));
pole_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(12, :))); 
% apex_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(13, :)));
% sliding_together = together_plot_filtered(:, logical(together_plot_filtered(14, :)));
scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond', 'MarkerEdgeColor','k'); hold on  % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'square', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', '^', 'MarkerEdgeColor','k'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', 'v', 'MarkerEdgeColor','k'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', 'pentagram', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(trapped_together(2, :)', trapped_together(5, :)', 200, trapped_together(1, :)', 'diamond','MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7]); hold on 
scatter(bypass_edge_together(2, :)', bypass_edge_together(5, :)', 200, bypass_edge_together(1, :)', 'Filled','MarkerEdgeColor','k'); hold on
scatter(bypass_tip_together(2, :)', bypass_tip_together(5, :)', 200, bypass_tip_together(1, :)', 'Filled', 'square','MarkerEdgeColor','k'); hold on
scatter(pole_vaulting_together(2, :)', pole_vaulting_together(5, :)', 200, pole_vaulting_together(1, :)', 'Filled', '^','MarkerEdgeColor','k'); hold on
% scatter(apex_vaulting_together(2, :)', apex_vaulting_together(5, :)', 200, apex_vaulting_together(1, :)', 'Filled', 'v','MarkerEdgeColor','k'); hold on
% scatter(sliding_together(2, :)', sliding_together(5, :)', 200, sliding_together(1, :)', 'Filled', 'pentagram', 'MarkerEdgeColor','k');
cmap(size(together_plot,2)); 
hcb=colorbar;
title(hcb,'$Deviation\ (\delta/h_{obs})$','FontSize', 20,'Interpreter', 'latex'); grid on
set(gca,'FontSize',16);
xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 24,'Interpreter', 'latex');
ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 24,'Interpreter', 'latex');
legend({'Trapping','Below','Above','Pole-vaulting'}, 'Location', 'southeast','FontSize', 14,'Interpreter', 'latex')
if range_L_low ~= 0; xlim([range_L_low range_L_up]); end
if range_y0_low ~= -10; ylim([range_y0_low range_y0_up]); end
% f=gcf;
% exportgraphics(f,'delta_vs_L-y0_classification_alldata_20230105.png','Resolution',100)


%% plot chi_0 & y_0 with classification and further data cleaning (optional)):  
figure('color', 'w'); set(gcf, 'Position', [100 100 1000 500]);
% cmap = cmocean('thermal');

together_plot_filtered = together_plot;

% together_plot_filtered(:, together_plot_filtered(10, :) > 800) = []; % remove the cases that are too fast
% together_plot_filtered(:, together_plot_filtered(10, :) < 400) = []; % remove the cases that are too slow
abs_delta_U = abs(together_plot_filtered(21, :) - together_plot_filtered(20, :)); % |U_f - u_0|
abs_delta_chi = abs(together_plot_filtered(16, :) - together_plot_filtered(3, :)); % |chi_f - chi_0|

together_plot_filtered(:, and(abs_delta_U > 100, abs_delta_chi < 10)) = []; % remove the cases that have large U-differences but very small chi-differences.

onlypass_ind = ~(logical(together_plot_filtered(12, :))); 
% the case index without pole-vaulting, apex-vaulting and sliding.
bypass_edge_together = together_plot_filtered(:, and(~logical(together_plot_filtered(8, :)), onlypass_ind));
bypass_tip_together = together_plot_filtered(:, and(logical(together_plot_filtered(8, :)), onlypass_ind));
pole_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(12, :))); 

plot(bypass_edge_together(3, :), bypass_edge_together(5, :), 'o','MarkerSize', 18,'MarkerEdgeColor','red','LineWidth',2); hold on
plot(bypass_tip_together(3, :), bypass_tip_together(5, :),  'square','MarkerSize', 20,'MarkerEdgeColor',[0 .5 0],'LineWidth',2); hold on
plot(pole_vaulting_together(3, :), pole_vaulting_together(5, :),  '^','MarkerSize', 18,'MarkerEdgeColor','blue','LineWidth',2); hold on
plot(trapped_together(3, :), trapped_together(5, :), 'diamond','MarkerSize', 18,'MarkerEdgeColor',[0.92, 0.70, 0.22],'LineWidth',2); hold on 

xlabel('$\theta_0$','FontSize', 18,'Interpreter', 'latex'); 
ylabel('$y_0/h_{\rm obs}$','FontSize', 18,'Interpreter', 'latex');

xlim([-10 10]); ylim([0 1]); 
set_plot(gcf, gca)

f=gcf;
savefig(f,'D:\Dropbox\Collaboration - LadHyX\Give_to_Clement\FSI - Rigid Fiber &  Individual Obstacle\delta_vs_theta0-y0_L0.5-1.5_20230629.fig')


%% plot the deviation vs chi_0 & y_0 (with classification and further data cleaning (optional)):  
figure('color', 'w'); set(gcf, 'Position', [100 100 1000 500]);
% cmap = cmocean('thermal');
cmap = colormap("jet");
together_plot_filtered = together_plot;
names_plot_filtered = names_plot;
% together_plot_filtered(:, together_plot_filtered(10, :) > 800) = []; % remove the cases that are too fast
% together_plot_filtered(:, together_plot_filtered(10, :) < 400) = []; % remove the cases that are too slow
abs_delta_U = abs(together_plot_filtered(21, :) - together_plot_filtered(20, :)); % |U_f - u_0|
abs_delta_chi = abs(together_plot_filtered(16, :) - together_plot_filtered(3, :)); % |chi_f - chi_0|

together_plot_filtered(:, and(abs_delta_U > 100, abs_delta_chi < 10)) = []; % remove the cases that have large U-differences but very small chi-differences.

onlypass_ind = ~(logical(together_plot_filtered(12, :))); 
% the case index without pole-vaulting, apex-vaulting and sliding.
bypass_edge_together = together_plot_filtered(:, and(~logical(together_plot_filtered(8, :)), onlypass_ind));
bypass_tip_together = together_plot_filtered(:, and(logical(together_plot_filtered(8, :)), onlypass_ind));
pole_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(12, :))); 
% apex_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(13, :)));
% sliding_together = together_plot_filtered(:, logical(together_plot_filtered(14, :)));
scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond', 'MarkerEdgeColor','k'); hold on  % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'square', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', '^', 'MarkerEdgeColor','k'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', 'v', 'MarkerEdgeColor','k'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', 'pentagram', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(trapped_together(3, :)', trapped_together(5, :)', 200, trapped_together(1, :)', 'diamond','MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7]); hold on 
scatter(bypass_edge_together(3, :)', bypass_edge_together(5, :)', 200, bypass_edge_together(1, :)', 'Filled','MarkerEdgeColor','k'); hold on
scatter(bypass_tip_together(3, :)', bypass_tip_together(5, :)', 200, bypass_tip_together(1, :)', 'Filled', 'square','MarkerEdgeColor','k'); hold on
scatter(pole_vaulting_together(3, :)', pole_vaulting_together(5, :)', 200, pole_vaulting_together(1, :)', 'Filled', '^','MarkerEdgeColor','k'); hold on
% scatter(apex_vaulting_together(3, :), apex_vaulting_together(5, :), 100, apex_vaulting_together(1, :), 'Filled', 'v','MarkerEdgeColor','k'); hold on
% scatter(sliding_together(3, :), sliding_together(5, :), 100, sliding_together(1, :), 'Filled', 'pentagram', 'MarkerEdgeColor','k');
cmap(size(together_plot,2)); 
hcb=colorbar;
title(hcb,'$Deviation\ (\delta/h_{obs})$','FontSize', 20,'Interpreter', 'latex'); grid on
set(gca,'FontSize',16);
xlabel('$Initial\ angle\ \theta_0\ (^{\circ})$','FontSize', 24,'Interpreter', 'latex');
ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 24,'Interpreter', 'latex');
% title('$0.45<L/l_{obs}<1.05$','FontSize', 20,'Interpreter', 'latex')
legend({'Trapping','Below','Above','Pole-vaulting'}, 'Location', 'southeast','FontSize', 14,'Interpreter', 'latex')
% xlim([-10 10]); ylim([0 1]); caxis([-0.3 0.3]);
% f=gcf;
% exportgraphics(f,'delta_vs_theta0-y0_L0.5-1.5_20230110.png','Resolution',100)



%% plot the y_0 vs L & deviation (with classification and further data cleaning (optional)):
figure('color', 'w'); set(gcf, 'Position', [100 100 1000 500]);
% cmap = cmocean('thermal');
cmap = colormap("jet");
together_plot_filtered = together_plot;
% together_plot_filtered(:, together_plot_filtered(10, :) > 800) = []; % remove the cases that are too fast
% together_plot_filtered(:, together_plot_filtered(10, :) < 400) = []; % remove the cases that are too slow
abs_delta_U = abs(together_plot_filtered(21, :) - together_plot_filtered(20, :)); % |U_f - u_0|
abs_delta_chi = abs(together_plot_filtered(16, :) - together_plot_filtered(3, :)); % |chi_f - chi_0|
together_plot_filtered(:, and(abs_delta_U > 100, abs_delta_chi < 10)) = []; % remove the cases that have large U-differences but very small chi-differences.
onlypass_ind = ~(logical(together_plot_filtered(12, :))); 
% the case index without pole-vaulting, apex-vaulting and sliding.
bypass_edge_together = together_plot_filtered(:, and(~logical(together_plot_filtered(8, :)), onlypass_ind));
bypass_tip_together = together_plot_filtered(:, and(logical(together_plot_filtered(8, :)), onlypass_ind));
pole_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(12, :))); 
% apex_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(13, :)));
% sliding_together = together_plot_filtered(:, logical(together_plot_filtered(14, :)));
% scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond', 'MarkerEdgeColor','k'); hold on  % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'square', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', '^', 'MarkerEdgeColor','k'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', 'v', 'MarkerEdgeColor','k'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', 'pentagram', 'MarkerEdgeColor','k'); hold on % for legend only
% scatter(trapped_together(2, :), trapped_together(1, :), 200, trapped_together(5, :), 'diamond','MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7]); hold on 
scatter(bypass_edge_together(2, :), bypass_edge_together(1, :), 200, bypass_edge_together(5, :), 'Filled','MarkerEdgeColor','k'); hold on
scatter(bypass_tip_together(2, :), bypass_tip_together(1, :), 200, bypass_tip_together(5, :), 'Filled', 'square','MarkerEdgeColor','k'); hold on
scatter(pole_vaulting_together(2, :), pole_vaulting_together(1, :), 200, pole_vaulting_together(5, :), 'Filled', '^','MarkerEdgeColor','k'); hold on
% scatter(apex_vaulting_together(2, :), apex_vaulting_together(1, :), 200, apex_vaulting_together(5, :), 'Filled', 'v','MarkerEdgeColor','k'); hold on
% scatter(sliding_together(2, :), sliding_together(1, :), 200, sliding_together(5, :), 'Filled', 'pentagram', 'MarkerEdgeColor','k');
cmap(size(together_plot,2)); 
hcb=colorbar;
title(hcb,'$Initial\ position\ (y_0/h_{obs})$','FontSize', 20,'Interpreter', 'latex'); grid on
set(gca,'FontSize',16);
xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 24,'Interpreter', 'latex');
ylabel('$Deviation\ (\delta/h_{obs})$','FontSize', 24,'Interpreter', 'latex');
legend({'Below','Above','Pole-vaulting'}, 'Location', 'southwest','FontSize', 14,'Interpreter', 'latex')
xlim([0.5 1.4]); ylim([-0.2 0.2]); caxis([0 1]);
% f=gcf;
% exportgraphics(f,'y0_vs_L-delta_classification_filtered_20221204.png','Resolution',100)
% [the_L, the_y0] = ginput(1); % pick up the point you want to show the trajectory.
% the_loc = intersect(find(together_plot(2, :)>the_L*0.98 & together_plot(2, :)<the_L*1.02),...
%     find(together_plot(1, :)>the_y0*1.02 & together_plot(1, :)<the_y0*0.98));  % Change the control range if there is an error.
% names_plot(the_loc)     % show the file

%% plot the y_0 vs deviation & (1-U_bar/U0) <interaction1> (with classification and further data cleaning (optional)):
figure('color', 'w'); set(gcf, 'Position', [100 100 1000 500]);
% cmap = cmocean('thermal');
cmap = colormap("jet");
together_plot_filtered = together_plot;
% together_plot_filtered(:, together_plot_filtered(10, :) > 800) = []; % remove the cases that are too fast
% together_plot_filtered(:, together_plot_filtered(10, :) < 400) = []; % remove the cases that are too slow
abs_delta_U = abs(together_plot_filtered(21, :) - together_plot_filtered(20, :)); % |U_f - u_0|
abs_delta_chi = abs(together_plot_filtered(16, :) - together_plot_filtered(3, :)); % |chi_f - chi_0|
together_plot_filtered(:, and(abs_delta_U > 100, abs_delta_chi < 10)) = []; % remove the cases that have large U-differences but very small chi-differences.
onlypass_ind = ~(logical(together_plot_filtered(12, :))); 
% the case index without pole-vaulting, apex-vaulting and sliding.
bypass_edge_together = together_plot_filtered(:, and(~logical(together_plot_filtered(8, :)), onlypass_ind));
bypass_tip_together = together_plot_filtered(:, and(logical(together_plot_filtered(8, :)), onlypass_ind));
pole_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(12, :))); 
apex_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(13, :)));
sliding_together = together_plot_filtered(:, logical(together_plot_filtered(14, :)));
scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond', 'MarkerEdgeColor','k'); hold on  % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'square', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', '^', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'v', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'pentagram', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(trapped_together(17, :), trapped_together(1, :), 200, trapped_together(5, :), 'diamond','MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7]); hold on 
scatter(bypass_edge_together(17, :), bypass_edge_together(1, :), 200, bypass_edge_together(5, :), 'Filled', 'MarkerEdgeColor','k'); hold on
scatter(bypass_tip_together(17, :), bypass_tip_together(1, :), 200, bypass_tip_together(5, :), 'Filled', 'square', 'MarkerEdgeColor','k'); hold on
scatter(pole_vaulting_together(17, :), pole_vaulting_together(1, :), 200, pole_vaulting_together(5, :), 'Filled', '^', 'MarkerEdgeColor','k'); hold on
scatter(apex_vaulting_together(17, :), apex_vaulting_together(1, :), 200, apex_vaulting_together(5, :), 'Filled', 'v', 'MarkerEdgeColor','k'); hold on
scatter(sliding_together(17, :), sliding_together(1, :), 200, sliding_together(5, :), 'Filled', 'pentagram', 'MarkerEdgeColor','k');
cmap(size(together_plot_filtered,2)); 
hcb=colorbar; caxis([0 1])
title(hcb,'$Initial\ position\ (y_0/h_{obs})$','FontSize', 20,'Interpreter', 'latex'); grid on
set(gca,'FontSize',16);
xlabel('$1-\bar{U}/U_0$','FontSize', 24,'Interpreter', 'latex');
ylabel('$Deviation\ (\delta/h_{obs})$','FontSize', 24,'Interpreter', 'latex');
legend({'Trapping','Passing Below','Passing Above (apex side)','Pole Vaulting','Apex Vaulting', ...
    'Sliding'}, 'Location', 'southwest','FontSize', 14,'Interpreter', 'latex')
% if range_L_low ~= 0; xlim([range_L_low range_L_up]); end
% if range_y0_low ~= -10; ylim([range_y0_low range_y0_up]); end
% f=gcf;
% exportgraphics(f,'y0_vs_interaction-delta_classification_datacleaning.png','Resolution',100)

%% plot the y_0 vs deviation & max(abs(U0-U(t)))/U0 <interaction2> (with classification and further data cleaning (optional)):
figure('color', 'w'); set(gcf, 'Position', [100 100 1000 500]);
% cmap = cmocean('thermal');
cmap = colormap("jet");
together_plot_filtered = together_plot;
% together_plot_filtered(:, together_plot_filtered(10, :) > 800) = []; % remove the cases that are too fast
% together_plot_filtered(:, together_plot_filtered(10, :) < 400) = []; % remove the cases that are too slow
abs_delta_U = abs(together_plot_filtered(21, :) - together_plot_filtered(20, :)); % |U_f - u_0|
abs_delta_chi = abs(together_plot_filtered(16, :) - together_plot_filtered(3, :)); % |chi_f - chi_0|
together_plot_filtered(:, and(abs_delta_U > 100, abs_delta_chi < 10)) = []; % remove the cases that have large U-differences but very small chi-differences.
onlypass_ind = ~(logical(together_plot_filtered(12, :))); 
% the case index without pole-vaulting, apex-vaulting and sliding.
bypass_edge_together = together_plot_filtered(:, and(~logical(together_plot_filtered(8, :)), onlypass_ind));
bypass_tip_together = together_plot_filtered(:, and(logical(together_plot_filtered(8, :)), onlypass_ind));
pole_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(12, :))); 
apex_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(13, :)));
sliding_together = together_plot_filtered(:, logical(together_plot_filtered(14, :)));
scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond', 'MarkerEdgeColor','k'); hold on  % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'square', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', '^', 'MarkerEdgeColor','k'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', 'v', 'MarkerEdgeColor','k'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', 'pentagram', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(trapped_together(18, :), trapped_together(1, :), 200, trapped_together(5, :), 'diamond','MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7]); hold on 
scatter(bypass_edge_together(18, :), bypass_edge_together(1, :), 200, bypass_edge_together(5, :), 'Filled', 'MarkerEdgeColor','k'); hold on
scatter(bypass_tip_together(18, :), bypass_tip_together(1, :), 200, bypass_tip_together(5, :), 'Filled', 'square', 'MarkerEdgeColor','k'); hold on
scatter(pole_vaulting_together(18, :), pole_vaulting_together(1, :), 200, pole_vaulting_together(5, :), 'Filled', '^', 'MarkerEdgeColor','k'); hold on
% scatter(apex_vaulting_together(18, :), apex_vaulting_together(1, :), 200, apex_vaulting_together(5, :), 'Filled', 'v', 'MarkerEdgeColor','k'); hold on
% scatter(sliding_together(18, :), sliding_together(1, :), 200, sliding_together(5, :), 'Filled', 'pentagram', 'MarkerEdgeColor','k');
cmap(size(together_plot_filtered,2)); 
hcb=colorbar; caxis([0 1])
title(hcb,'$Initial\ position\ (y_0/h_{obs})$','FontSize', 20,'Interpreter', 'latex'); grid on
set(gca,'FontSize',16);
xlabel('$max\left|U_0-U(t)\right|/U_0$','FontSize', 24,'Interpreter', 'latex');
ylabel('$Deviation\ (\delta/h_{obs})$','FontSize', 24,'Interpreter', 'latex');
legend({'Trapping','Below','Above','Pole-vaulting'...
    }, 'Location', 'southwest','FontSize', 14,'Interpreter', 'latex')
% title('$-10<{\theta}_0<10,\ 0.5<L<1$','FontSize', 22,'Interpreter', 'latex')
% if range_L_low ~= 0; xlim([range_L_low range_L_up]); end
% if range_y0_low ~= -10; ylim([range_y0_low range_y0_up]); end
% f=gcf;
% exportgraphics(f,'y0_vs_interaction2-delta_classification_datacleaning_Chi0-m10to10_L-0p5to1.png','Resolution',100)
[the_inter2, the_delta] = ginput(1); % pick up the point you want to show the trajectory.
the_loc = intersect(find(together_plot(18, :)>=the_inter2*0.98 & together_plot(18, :)<=the_inter2*1.02),...
    find(together_plot(1, :)>=the_delta*1.02 & together_plot(1, :)<=the_delta*0.98));  % Change the control range if there is an error.
names_plot(the_loc)     % show the file

%% plot the y_0 vs deviation & <interaction3> (with classification and further data cleaning (optional)):
figure('color', 'w'); set(gcf, 'Position', [100 100 1000 500]);
% cmap = cmocean('thermal');
cmap = colormap("jet");
together_plot_filtered = together_plot;
% together_plot_filtered(:, together_plot_filtered(10, :) > 800) = []; % remove the cases that are too fast
% together_plot_filtered(:, together_plot_filtered(10, :) < 400) = []; % remove the cases that are too slow
abs_delta_U = abs(together_plot_filtered(21, :) - together_plot_filtered(20, :)); % |U_f - u_0|
abs_delta_chi = abs(together_plot_filtered(16, :) - together_plot_filtered(3, :)); % |chi_f - chi_0|
together_plot_filtered(:, and(abs_delta_U > 100, abs_delta_chi < 10)) = []; % remove the cases that have large U-differences but very small chi-differences.
onlypass_ind = ~(logical(together_plot_filtered(12, :))); 
% the case index without pole-vaulting, apex-vaulting and sliding.
bypass_edge_together = together_plot_filtered(:, and(~logical(together_plot_filtered(8, :)), onlypass_ind));
bypass_tip_together = together_plot_filtered(:, and(logical(together_plot_filtered(8, :)), onlypass_ind));
pole_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(12, :))); 
apex_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(13, :)));
sliding_together = together_plot_filtered(:, logical(together_plot_filtered(14, :)));
scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond', 'MarkerEdgeColor','k'); hold on  % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'square', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', '^', 'MarkerEdgeColor','k'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', 'v', 'MarkerEdgeColor','k'); hold on % for legend only
% scatter(nan, nan, 1, nan, 'filled', 'k', 'pentagram', 'MarkerEdgeColor','k'); hold on % for legend only
scatter(trapped_together(19, :), trapped_together(1, :), 200, trapped_together(5, :), 'diamond','MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7]); hold on 
scatter(bypass_edge_together(19, :), bypass_edge_together(1, :), 200, bypass_edge_together(5, :), 'Filled','MarkerEdgeColor','k'); hold on
scatter(bypass_tip_together(19, :), bypass_tip_together(1, :), 200, bypass_tip_together(5, :), 'Filled', 'square', 'MarkerEdgeColor','k'); hold on
scatter(pole_vaulting_together(19, :), pole_vaulting_together(1, :), 200, pole_vaulting_together(5, :), 'Filled', '^', 'MarkerEdgeColor','k'); hold on
% scatter(apex_vaulting_together(19, :), apex_vaulting_together(1, :), 200, apex_vaulting_together(5, :), 'Filled', 'v', 'MarkerEdgeColor','k'); hold on
% scatter(sliding_together(19, :), sliding_together(1, :), 200, sliding_together(5, :), 'Filled', 'pentagram', 'MarkerEdgeColor','k');
cmap(size(together_plot_filtered,2)); 
hcb=colorbar; caxis([0 1])
title(hcb,'$Initial\ position\ (y_0/h_{obs})$','FontSize', 20,'Interpreter', 'latex'); grid on
set(gca,'FontSize',16);
xlabel('$Contact\ probability$','FontSize', 24,'Interpreter', 'latex');  % for interaction3
ylabel('$Deviation\ (\delta/h_{obs})$','FontSize', 24,'Interpreter', 'latex');
legend({'Trapping','Below','Above','Pole-vaulting'...
    }, 'Location', 'southwest','FontSize', 14,'Interpreter', 'latex')
% title('$0.5<L<1$','FontSize', 22,'Interpreter', 'latex')
% if range_L_low ~= 0; xlim([range_L_low range_L_up]); end
% if range_y0_low ~= -10; ylim([range_y0_low range_y0_up]); end
% f=gcf;
% exportgraphics(f,'y0_vs_interaction3-delta_classification_datacleaning_L-0p5to1.png','Resolution',100)
% [the_inter2, the_delta] = ginput(1); % pick up the point you want to show the trajectory.
% the_loc = intersect(find(together_plot(19, :)>=the_inter2*0.98 & together_plot(19, :)<=the_inter2*1.02),...
%     find(together_plot(1, :)>=the_delta*1.02 & together_plot(1, :)<=the_delta*0.98));  % Change the control range if there is an error.
% names_plot(the_loc)     % show the file

%% plot the speed vs L & y_0: 
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% cmap = cmocean('thermal');
cmap = colormap("jet");
scatter(together_plot(2, :), together_plot(5, :), 100, together_plot(10, :), 'Filled', 'hexagram', 'MarkerEdgeColor','k') 
cmap(size(together_plot,2)); 
hcb=colorbar;
title(hcb,'$Average\ speed\ (\mu{m}/s)$','FontSize', 16,'Interpreter', 'latex');
grid on
set(gca,'FontSize',16);
xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,'Ubar_vs_L-y0.png','Resolution',100)

%% plot (1-U_bar/U0) <interaction1> vs L & y_0:
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% cmap = cmocean('thermal');
cmap = colormap("jet");
scatter(together_plot(2, :), together_plot(5, :), 100, together_plot(17, :), 'Filled', 'hexagram', 'MarkerEdgeColor','k')
cmap(size(together_plot,2)); 
hcb=colorbar;
title(hcb,'$1-\bar{U}/U_0$','FontSize', 16,'Interpreter', 'latex');
grid on
set(gca,'FontSize',16);
xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,'interaction_vs_L-y0.png','Resolution',100)

%% plot the delta Chi vs L & y_0:
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% cmap = cmocean('thermal');
cmap = colormap("jet");
scatter(together_plot(2, :), together_plot(5, :), 100, together_plot(11, :)', 'Filled', 'hexagram', 'MarkerEdgeColor','k')    % add ' to avoid obfuscating to a triplet.
cmap(size(together_plot,2)); 
hcb=colorbar;
title(hcb,'$\theta_f - \theta_0\ (^{\circ})$','FontSize', 16,'Interpreter', 'latex');
grid on
set(gca,'FontSize',16);
xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,'deltaChi_vs_L-y0.png','Resolution',100)



%%%%%%%%%%%%%%%%%%%%%%%%% other plots %%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the |delta Chi| vs deviation 
figure('color', 'w'); set(gcf, 'Position', [100 100 600 400]);
plot(abs(together_plot(11, :)), together_plot(1, :), 'Color','r', 'LineStyle','none', 'Marker','*', 'MarkerSize', 5)
set(gca,'FontSize',16);
xlabel('$\left|\theta_f - \theta_0 \right| (^{\circ})$','FontSize',22,'Interpreter', 'latex');
ylabel('$Deviation\ (\delta/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,'absDeltaChi_vs_delta.png','Resolution',100)

%% plot the (1 - U_bar/U0) <interaction1>  vs deviation
figure('color', 'w'); set(gcf, 'Position', [100 100 600 400]);
plot(together_plot(17, :), together_plot(1, :), 'Color','r', 'LineStyle','none', 'Marker','*', 'MarkerSize', 5)
set(gca,'FontSize',16);
xlabel('$1-\bar{U}/U_0$','FontSize',22,'Interpreter', 'latex');
ylabel('$Deviation\ (\delta/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,'interaction_vs_delta.png','Resolution',100)

%% plot the (speed_downstream-speed_upstream) vs Chi_0/Chi_f
figure('color', 'w'); set(gcf, 'Position', [100 100 600 400]);
yyaxis left
plot((together_plot(21, :) - together_plot(20, :)), together_plot(3, :), 'LineStyle','none', 'Marker','*', 'MarkerSize', 7)
ylabel('$\theta_0\ (^{\circ})$','FontSize', 22,'Interpreter', 'latex');
yyaxis right
plot((together_plot(21, :) - together_plot(20, :)), together_plot(16, :), 'LineStyle','none', 'Marker','o', 'MarkerSize', 7)
ylabel('$\theta_f\ (^{\circ})$','FontSize', 22,'Interpreter', 'latex');
set(gca,'FontSize',16);
xlabel('$U_{downstream}-U_{upstream}\ (\mu{m}/s)$','FontSize',22,'Interpreter', 'latex');
xlim([-300 300])
% f=gcf;
% exportgraphics(f,'deltaU_vs_chi0__eltaU_vs_chif.png','Resolution',100)

%% plot <interaction3> 'contact probability' vs. deviation (with further data cleaning):
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
together_plot_filtered = together_plot;
abs_delta_U = abs(together_plot_filtered(21, :) - together_plot_filtered(20, :)); % |U_f - u_0|
abs_delta_chi = abs(together_plot_filtered(16, :) - together_plot_filtered(3, :)); % |chi_f - chi_0|
together_plot_filtered(:, and(abs_delta_U > 100, abs_delta_chi < 10)) = []; % remove the cases that have large U-differences but very small chi-differences.

plot(together_plot_filtered(19, :), together_plot_filtered(1, :), 'ok', 'MarkerSize', 9, 'LineWidth', 1.5)

xlim([-0.1 1.1]); ylim([-0.2 0.3]);
xlabel('Contact probability','FontSize', 24,'FontName', 'Times New Roman'); 
ylabel('$\delta$','FontSize', 24,'Interpreter', 'latex');
text(-0.05, 0.25, 'Experiment','FontSize', 24, 'Interpreter', 'latex','BackgroundColor',[.7 .7 .7])
set(gca,'Box', 'On','XGrid', 'On','YGrid', 'On','FontSize', 24,'TickLabelInterpreter','latex')

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['F:\Processing & Results\FSI - Rigid Fiber &  Individual ' ...
    'Obstacle\Figures\about interaction index\Contact_probability_delta_exp.pdf']);



%%
function data_out = delta_correction(data_in, slope)
% Only for obstacle points below!
% data_in should be 'All_data'.
% slope is a number (notice about the sign!).
data_out = data_in;
for foo = 1:length(data_in.delta_y)
    data_out.delta_y(foo) = data_in.delta_y(foo) - slope * (data_in.final_x(foo) - data_in.initial_x(foo));
end
end