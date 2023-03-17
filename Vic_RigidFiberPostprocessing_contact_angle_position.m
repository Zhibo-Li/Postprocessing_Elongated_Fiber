clear; close all; clc;

% basic information.
h_obs = 75; l_obs = 86.6; Obj_Mag = 0.63;

% laod data and put all in one.
load('F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\20220624-SU8_Fibers-Individual_triangularPillar_uppoint_contact_information.mat');
All_data1 = delta_correction(All_data, 0);
load('F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\20220913-SU8_Fibers-Individual_triangularPillar_uppoint_contact_information.mat');
All_data2 = All_data;
load('F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\20221004-SU8_Fibers-Individual_triangularPillar_uppoint_contact_information.mat');
All_data3 = All_data;
load('F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\20221005-SU8_Fibers-Individual_triangularPillar_uppoint_contact_information.mat');
All_data4 = All_data;
load('F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\20230102-SU8_Fibers-Individual_triangularPillar_uppoint_contact_information.mat');
All_data5 = delta_correction(All_data, 0);
All_data = [All_data1, All_data2, All_data3, All_data4, All_data5];

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
Chi_c = [All_data(:).Chi_c]; % interaction index 1
y_c = [All_data(:).y_c]; % interaction index 2

% calculate the delta chi (between chi_f and chi_0)
Chi = [All_data(:).Chi];
for i = 1:length(Chi)
    dt = diff(timestamps{1, i}); % [time(i+1) - time(i)]
    dchi = diff(Chi{1, i})'; % chi(i+1) - chi(i)
%     dchi_dt = dchi ./ (dt * 100); % [chi(i+1) - chi(i)] / [time(i+1) - time(i)]. Here * 100 for convenience.
    n_plus = length(find(dchi<-60));  n_minus = length(find(dchi>60)); % check if there is a rotation and determine the direction via the value of d_chi
    delta_Chi(i) = Chi{1, i}(end) - Chi{1, i}(1) + n_plus*180 - n_minus*180;
    Chi_f(i) = - Chi{1, i}(end); % the orientation of the fiber in the last frame (most downstream)
end
delta_Chi(18) = delta_Chi(18) + 180; % correct these values.
delta_Chi(87) = delta_Chi(87)  + 180;
delta_Chi(102) = delta_Chi(102) + 180;
delta_Chi(110) = delta_Chi(110) + 180;
delta_Chi = -delta_Chi; % delta_chi: counterclockwise is positive (meet the definition).

% put all the information in a matrix
together = [norm_delta_y; norm_contourL; Chi_0; norm_initial_x; norm_initial_y;...
    norm_final_x; norm_final_y; if_bypass_tip_together; ave_speed; speed_ave_all; ...
    delta_Chi; PoleVaulting; ApexVaulting; Sliding; Trapping; ... 
    Chi_f; Chi_c; y_c; speed_upstream; speed_downstream];
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
% %     No.17 row: chi_c
% %     No.18 row: y_c
% %     No.19 row: upstream speed
% %     No.20 row: downstream speed

names(~logical(acceptability)) = [];
together(:, ~logical(acceptability)) = [];
Trapping(~logical(acceptability)) = [];
% these cases are not in the channel mid-plane by manual selection.
% In the excel of PoleVaultingCheck.xls, their 'Deviation correction based on tracers' are NaN.

trapped_names = names(logical(Trapping));  % extract the trapping case names
trapped_together = together(:, logical(Trapping)); % extract the trapping case information
names(logical(Trapping)) = [];
together(:, logical(Trapping)) = []; % remove trapping cases

% select reasonable cases.
range_upstream = -3; range_downstream = 3; % set acceptable range (norm_initial_x and norm_final_x)
names(together(4, :) > range_upstream) = [];
together(:, together(4, :) > range_upstream) = []; % remove the cases too close to the obstacle on the upstream side.
names(together(6, :) < range_downstream) = [];
together(:, together(6, :) < range_downstream) = []; % remove the cases too close to the obstacle on the downstream side.

% create dialog box to gather user input for plotting;
names_plot = names; together_plot = together;
prompt = {'The lower bound of the initial angle:', 'The upper bound of the initial angle:', ...
    'The lower bound of the contour length:','The upper bound of the contour length:'...
    'The lower bound of the initial position:','The upper bound of the initial position:'};
definput = {'-10', '10', 'nan', 'nan', '0', '1'};
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


% plot the chi_c vs  y_c (with classification and further data cleaning (optional)):
figure('color', 'w'); set(gcf, 'Position', [100 100 420 315]);
% cmap = cmocean('thermal');
cmap = colormap("jet");
together_plot_filtered = together_plot;
% together_plot_filtered(:, together_plot_filtered(10, :) > 800) = []; % remove the cases that are too fast
% together_plot_filtered(:, together_plot_filtered(10, :) < 400) = []; % remove the cases that are too slow
abs_delta_U = abs(together_plot_filtered(20, :) - together_plot_filtered(19, :)); % |U_f - u_0|
abs_delta_chi = abs(together_plot_filtered(16, :) - together_plot_filtered(3, :)); % |chi_f - chi_0|
together_plot_filtered(:, and(abs_delta_U > 100, abs_delta_chi < 10)) = []; % remove the cases that have large U-differences but very small chi-differences.
onlypass_ind = ~(logical(together_plot_filtered(12, :))); 
% the case index without pole-vaulting, apex-vaulting and sliding.
bypass_edge_together = together_plot_filtered(:, and(~logical(together_plot_filtered(8, :)), onlypass_ind));
bypass_tip_together = together_plot_filtered(:, and(logical(together_plot_filtered(8, :)), onlypass_ind));
pole_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(12, :))); 
% apex_vaulting_together = together_plot_filtered(:, logical(together_plot_filtered(13, :)));
% sliding_together = together_plot_filtered(:, logical(together_plot_filtered(14, :)));

% for legend
plot(nan, nan, 'diamond','MarkerSize', 8,'MarkerEdgeColor','k','MarkerFaceColor',[0.92, 0.70, 0.22]); 
plot(nan, nan, 'o','MarkerSize', 8,'MarkerEdgeColor','k','MarkerFaceColor','red');
plot(nan, nan,  'square','MarkerSize', 8,'MarkerEdgeColor','k','MarkerFaceColor', [0 .5 0]); 
plot(nan, nan,  '^','MarkerSize', 8,'MarkerEdgeColor','k','MarkerFaceColor','blue');  
% for plot
plot(trapped_together(17, :), trapped_together(18, :), 'diamond','MarkerSize', 10,'MarkerEdgeColor',[0.92, 0.70, 0.22],'LineWidth',2); hold on 
plot(bypass_edge_together(17, :), bypass_edge_together(18, :), 'o','MarkerSize', 10,'MarkerEdgeColor','red','LineWidth',2); hold on
plot(bypass_tip_together(17, :), bypass_tip_together(18, :),  'square','MarkerSize', 10,'MarkerEdgeColor',[0 .5 0],'LineWidth',2); hold on
plot(pole_vaulting_together(17, :), pole_vaulting_together(18, :),  '^','MarkerSize', 10,'MarkerEdgeColor','blue','LineWidth',2); hold on

xlabel('$\theta_c$','FontSize', 18,'Interpreter', 'latex'); 
ylabel('$y_c/h_{obs}$','FontSize', 18,'Interpreter', 'latex');
% title_txt = ['$-90 < \theta_0 < 90$'];
% title(title_txt,'FontSize', 18,'Interpreter', 'latex');
xlim([-50 50]); ylim([0 1]);
xticks([-45 -30 -15 0 15 30 45])
% legend({'Trapping','Below','Above','Pole-vaulting'}, 'Location', 'northwest','FontSize', 14,'Interpreter', 'latex')

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\Paper\contact information vs initial condition\thetaC-yC_dyn_exp_theta_0m10to10.png'],'Resolution',100)

% [the_L, the_y0] = ginput(1); % pick up the point you want to show the trajectory.
% the_loc = intersect(find(together_plot(17, :)<the_L*0.98 & together_plot(17, :)>the_L*1.02),...
%     find(together_plot(18, :)>the_y0*0.98 & together_plot(18, :)<the_y0*1.02));  % Change the control range if there is an error.
% names_plot(the_loc)     % show the file




%%%%%%%%%%% functions
function data_out = delta_correction(data_in, slope)
% Only for obstacle points below!
% data_in should be 'All_data'.
% slope is a number (notice about the sign!).
data_out = data_in;
for foo = 1:length(data_in.delta_y)
    data_out.delta_y(foo) = data_in.delta_y(foo) - slope * (data_in.final_x(foo) - data_in.initial_x(foo));
end
end