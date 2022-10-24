clear; close all; clc;
pathname = uigetdir('F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle', 'Choose a folder');  % input file
Files = dir(fullfile(pathname, 'results with timestamps', '*.mat'));
bg = imread(fullfile(pathname, 'results', "The Pillar.tif"));
the_excel = readcell(fullfile(pathname, 'PoleVaultingCheck.xls'), 'NumHeaderLines', 1);
avi_names = the_excel(:, 1);
im_BW = imbinarize(bg);
LL = regionprops(im_BW, 'Centroid');
Pillar_CoM  = LL.Centroid;
the_pillar_im = imcomplement(imfill(im_BW, 'holes'));
% cmap = cmocean('thermal');
Obj_Mag = 0.63;

im_no = 1:length(Files);
for ii = 1:length(Files)

    load([Files(ii).folder, filesep, Files(ii).name]);
    avi_No = find(cellfun(@(x) any(strcmp(x, append(extractBetween(Files(ii).name, ...
        'trajectory_', '_AABGR'), '.avi'))), avi_names));
    All_data.PoleVaulting(ii) = the_excel{avi_No, 2};
    All_data.Trapping(ii) = the_excel{avi_No, 3};
    All_data.Sliding(ii) = the_excel{avi_No, 4};
    All_data.ApexVaulting(ii) = the_excel{avi_No, 5};

%     %%% Plot trajectories
%     the_time = Good_case_frm_time * 100;
%     the_time = round(the_time - min(the_time));
%     
%     figure('color', 'w'); set(gcf, 'Position', [300 300 600 100]);
%     imshow(the_pillar_im(801:1200, :)); hold on
%     %%% define the range of loop indicator k based on the colormap
% %     if max(the_time) > 255
% %         k_range = find(the_time > 255, 1) - 1;
% %         for k =1:min(k_range, length(Good_case_frm))
% %             plot(xy.spl{Good_case_frm(k)}(:,1),2048-xy.spl{Good_case_frm(k)}(:,2)-800,'Color',cmap(the_time(k)+1,:),'LineWidth',1.5)
% %             hold on
% %         end
% %     else
% %         for k =1:length(Good_case_frm)
% %             plot(xy.spl{Good_case_frm(k)}(:,1),2048-xy.spl{Good_case_frm(k)}(:,2)-800,'Color',cmap(the_time(k)+1,:),'LineWidth',1.5)
% %             hold on
% %         end
% %     end
%     %%% draw trajectories with looping the colormap
%     for k =1:length(Good_case_frm)
%         plot(xy.spl{Good_case_frm(k)}(:,1),2048-xy.spl{Good_case_frm(k)}(:,2)-800,'Color',cmap(mod(the_time(k)+1, 256)+1,:),'LineWidth',1.5)
%         hold on
%     end
%     axis equal
%     hold off
% %     hcb=colorbar;
% %     title(hcb,'$Time\ (s)$','FontSize', 16,'Interpreter', 'latex');
%     f=gcf;
%     exportgraphics(f,[pathname,'\results with timestamps\',Files(ii).name,'_looping_color.png'],'Resolution',100)
%     close

    % add the names
    All_data.filename{ii} = Files(ii).name;    

    % add the timestamps
    All_data.timestamps{ii} = Good_case_frm_time - min(Good_case_frm_time);

    % add the contour lengths 
    sorted_lengths = sort(xy.arclen_spl(Good_case_frm));
    All_data.contour_length(ii) = mean(sorted_lengths) * Obj_Mag;  % Average to get the contour length (UNIT: um).

    % add the position information
    centroidxy = reshape(cell2mat(xy.centroid),2,numel(xy.centroid));
    centroidxy = centroidxy(:, Good_case_frm);
    All_data.delta_y(ii) = (centroidxy(2, 1) - centroidxy(2,end)) * Obj_Mag;
    All_data.initial_y(ii) = ((2048-centroidxy(2, 1)) - Pillar_CoM(2)) * Obj_Mag + 25; % Definition of y0 is the distance between CoM and the edge, so +25.
    All_data.initial_x(ii) = (centroidxy(1, 1) - Pillar_CoM(1)) * Obj_Mag;
    All_data.final_y(ii) = ((2048-centroidxy(2, end)) - Pillar_CoM(2)) * Obj_Mag + 25; % Definition of y0 is the distance between CoM and the edge, so +25.
    All_data.final_x(ii) = (centroidxy(1, end) - Pillar_CoM(1)) * Obj_Mag;
    All_data.CoM{ii} = centroidxy;
    %     hold on; plot(centroidxy(1,:), 2048-centroidxy(2,:)-800, "Color",'k', 'LineStyle','--');

    % calculate the average speed along x-direction (UNIT: um/s)
    All_data.ave_speed(ii) = ((centroidxy(1, end) - Pillar_CoM(1)) - (centroidxy(1, 1) - Pillar_CoM(1))) * Obj_Mag / (max(Good_case_frm_time) - min(Good_case_frm_time));

    % check if the fiber bypass the obtacle tip
    if 2048 - mean(centroidxy(2, and(Pillar_CoM(1)-100 < centroidxy(1, :), centroidxy(1, :) < Pillar_CoM(1)+100))) > Pillar_CoM(2) % if it goes below (in image system).
        bypass_tip = true;
    else
        bypass_tip = false;
    end
    All_data.bypass_tip(ii) = bypass_tip;

    % calculate the initial angle Chi_0
    crd = xy.crd{1,Good_case_frm(1)};
    centroidxy_k = xy.centroid{1,Good_case_frm(1)};
    Gyr = 1/size(crd,1) * [sum((crd(:, 1)-centroidxy_k(1)).^2),  sum((crd(:, 1)-centroidxy_k(1)) .* (crd(:, 2)-centroidxy_k(2)));
        sum((crd(:, 2)-centroidxy_k(2)) .* (crd(:, 1)-centroidxy_k(1))), sum((crd(:, 2)-centroidxy_k(2)).^2)];

    [eigenV,eigenD] = eig(Gyr);
    [d,ind] = sort(diag(eigenD));
    Ds = eigenD(ind,ind);
    Vs = eigenV(:,ind);
    All_data.Chi_0(ii)  = atan(Vs(2,2)/Vs(1,2))/pi*180;

    % calculate all the angles
    for jj = 1:length(Good_case_frm)
        crd = xy.crd{1,Good_case_frm(jj)};
        centroidxy_k = xy.centroid{1,Good_case_frm(jj)};
        Gyr = 1/size(crd,1) * [sum((crd(:, 1)-centroidxy_k(1)).^2),  sum((crd(:, 1)-centroidxy_k(1)) .* (crd(:, 2)-centroidxy_k(2)));
            sum((crd(:, 2)-centroidxy_k(2)) .* (crd(:, 1)-centroidxy_k(1))), sum((crd(:, 2)-centroidxy_k(2)).^2)];

        [eigenV,eigenD] = eig(Gyr);
        [d,ind] = sort(diag(eigenD));
        Ds = eigenD(ind,ind);
        Vs = eigenV(:,ind);
        Chi(jj) = atan(Vs(2,2)/Vs(1,2))/pi*180;
    end
    All_data.Chi{ii}  = Chi;
    clearvars Chi

end

% save([pathname, '_infos_contourL-fromAVG_with_all-Chi_timestamps.mat'], 'All_data');



%% Ploting
clear; close all; clc;

h_obs = 75; l_obs = 86.6; Obj_Mag = 0.63;

load('20220624-SU8_Fibers-Individual_triangularPillar_uppoint_infos_contourL-fromAVG_with_all-Chi_timestamps.mat');
All_data1 = All_data;
load('20220913-SU8_Fibers-Individual_triangularPillar_uppoint_infos_contourL-fromAVG_with_all-Chi_timestamps.mat');
All_data2 = All_data;
load('20221004-SU8_Fibers-Individual_triangularPillar_uppoint_infos_contourL-fromAVG_with_all-Chi_timestamps.mat');
All_data3 = All_data;
load('20221005-SU8_Fibers-Individual_triangularPillar_uppoint_infos_contourL-fromAVG_with_all-Chi_timestamps.mat');
All_data4 = All_data;
All_data = [All_data1, All_data2, All_data3, All_data4];

names = [All_data(:).filename];

delta_y = [All_data(:).delta_y];
norm_delta_y = delta_y / h_obs; % normalized delta
contourL = [All_data(:).contour_length];
norm_contourL = contourL / l_obs; % normalized L
Chi_0 = [All_data(:).Chi_0]; 

initial_x = [All_data(:).initial_x];
norm_initial_x = initial_x / h_obs; % normalized x_0
initial_y = [All_data(:).initial_y];
norm_initial_y = initial_y / h_obs; % normalized y_0
final_x = [All_data(:).final_x];
norm_final_x = final_x / h_obs; % normalized x_f
final_y = [All_data(:).final_y];
norm_final_y = final_y / h_obs; % normalized y_f
bypass_tip = [All_data(:).bypass_tip];
speed = [All_data(:).ave_speed];
timestamps = [All_data(:).timestamps];

%%% calculate the delta chi
Chi = [All_data(:).Chi];
% plot the chi evolution.
for i = 1:length(Chi)
    dt = diff(timestamps{1, i}); % [time(i+1) - time(i)]
    dchi = diff(Chi{1, i})'; % chi(i+1) - chi(i)
    dchi_dt = dchi ./ (dt * 100); % [chi(i+1) - chi(i)] / [time(i+1) - time(i)]. Here * 100 for convenience.

%     figure('color', 'w'); set(gcf, 'Position', [-2000 500 2000 300]);
%     subplot(1,4,1); plot(timestamps{1, i}, Chi{1, i}); ylim([-90 90]); title('chi') % plot chi vs. time.
%     subplot(1,4,2); plot(dt); ylim([0 0.1]); title('d_t')
%     subplot(1,4,3); plot(dchi); ylim([-180 180]); title('d_c_h_i')
%     hold on; line([0 300],[60 60],'color','b');
%     hold on; line([0 300],[-60 -60],'color','b'); hold off 
%     subplot(1,4,4); plot(dchi_dt); ylim([-15 15]); title('d_c_h_i/d_t')

    n_plus = length(find(dchi<-60));  n_minus = length(find(dchi>60)); % check if there is a rotation and determine the direction via the value of d_chi

%     if n_plus == 0 && n_minus == 0 && max(abs(dchi)) > 60
%         n_plus = n_plus + sum(dt(dchi<-60) > 0.1);  
%         n_minus = n_minus + sum(dt(dchi>60) > 0.1); 
%     end

    % delta_chi: counterclockwise is positive.
    delta_Chi(i) = Chi{1, i}(end) - Chi{1, i}(1) + n_plus*180 - n_minus*180;
end
delta_Chi(18) = delta_Chi(18) + 180; % correct these values.
delta_Chi(27) = delta_Chi(27) + 180;
delta_Chi(38) = delta_Chi(38) - 180;
delta_Chi(39) = delta_Chi(39) + 180;
delta_Chi(87) = delta_Chi(87)  + 180;
delta_Chi(102) = delta_Chi(102) + 180;
delta_Chi(104) = delta_Chi(104) + 180;
delta_Chi(110) = delta_Chi(110) + 180;

%%% calculate the dx/dt
CoM = [All_data(:).CoM];
for j = 1:length(CoM)
    dt = diff(timestamps{1, j}); % [time(i+1) - time(i)]
    dx = diff(CoM{1, j}(1, :))'; % [x(i+1) - x(i)]
    dx_dt = movmean(dx, 7) ./ movmean(dt, 7) * Obj_Mag; % [chi(i+1) - chi(i)] / [time(i+1) - time(i)]. Here * 100 for convenience.

%     figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
%     plot(timestamps{1, j}(1:end-1), dx_dt, 'Color','m', 'LineStyle','none', 'Marker','.', 'MarkerSize', 20); 
%     xlim([0 3]); ylim([0 1200])
%     xlabel('$Time\ (s)$','FontSize', 22,'Interpreter', 'latex');
%     ylabel('$U_x\ (\mu{m}/s)$','FontSize', 22,'Interpreter', 'latex');
%     f=gcf;
%     exportgraphics(f,[num2str(j),'_',names{1, j}(1:end-4)  ,'_U_x.png'],'Resolution',100)
%     close all

    num_upper_obs = sum(CoM{1, j}(1, :) < 300);
    speed_upstream(j) = mean(dx_dt(1:num_upper_obs));
    speed_ave_all(j) = mean(dx_dt);
end



% %% statistics
% figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% edges = 20:20:180;
% histogram(contourL,edges);
% set(gca,'FontSize',16);
% xlabel('$Contour\ length\ L\ ({\mu}m)$','FontSize', 22,'Interpreter', 'latex');
% ylabel('$Number\ of\ cases$','FontSize', 22,'Interpreter', 'latex');
% xlim([20 160]);
% % f=gcf;
% % exportgraphics(f,'Statistics_contourL.png','Resolution',100)
% 
% %% statistics
% figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% edges = -90:20:90;
% histogram(Chi_0,edges);
% set(gca,'FontSize',16);
% xlabel('$Initial\ angle\ \chi_0\ (^{\circ})$','FontSize', 22,'Interpreter', 'latex');
% ylabel('$Number\ of\ cases$','FontSize', 22,'Interpreter', 'latex');
% xlim([-90 90]);
% % f=gcf;
% % exportgraphics(f,'Statistics_Chi0.png','Resolution',100)
% 
% %% statistics
% figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% edges = -0.2:0.2:1.2;
% histogram(norm_initial_y,edges);
% set(gca,'FontSize',16);
% xlabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% ylabel('$Number\ of\ cases$','FontSize', 22,'Interpreter', 'latex');
% xlim([-0.2 1.2]);
% % f=gcf;
% % exportgraphics(f,'Statistics_y0.png','Resolution',100)

together = [norm_delta_y; norm_contourL; Chi_0; norm_initial_x; norm_initial_y; norm_final_x; norm_final_y; bypass_tip; ones(1, length(CoM))-speed_ave_all./speed_upstream; delta_Chi];
% No.9 row: ones(1, length(CoM))-speed./speed_upstream: 1 - U_bar/U0.
trapped_names = names(together(6, :) < 0);  % extract the trapping case names
trapped_together = together(:, together(6, :) < 0); % extract the trapping case information
names(together(4, :) > -7) = [];
together(:, together(4, :) > -7) = []; % remove the cases too close to the obstacle on the upstream side.
names(together(6, :) < 7) = [];
together(:, together(6, :) < 7) = []; % remove the cases too close to the obstacle on the downstream side.

%%% the map: without any 'filter'.
% figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% cmap = cmocean('thermal');
% scatter(together(2, :), together(5, :), 50, together(1, :), 'Filled')
% cmap(size(contourL,2)); 
% hcb=colorbar;
% title(hcb,'$Deviation\ (\delta/h_{obs})$','FontSize', 16,'Interpreter', 'latex');
% grid on
% set(gca,'FontSize',16);
% xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
% ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% % f=gcf;
% % exportgraphics(f,'The_map_full_AVG-L.png','Resolution',100)

range_chi_low = -10; range_chi_up = 10;
names(together(3, :) < range_chi_low) = [];
together(:, together(3, :) < range_chi_low) = [];  % the initial angle: -10 < Chi_0 < 10.
names(together(3, :) > range_chi_up) = [];
together(:, together(3, :) > range_chi_up) = [];
%%% the map -10 < Chi_0 < 10.
% figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% cmap = cmocean('thermal');
% scatter(together(2, :), together(5, :), 50, together(1, :), 'Filled')
% cmap(size(contourL,2)); 
% hcb=colorbar;
% title(hcb,'$Deviation\ (\delta/h_{obs})$','FontSize', 16,'Interpreter', 'latex');
% grid on
% set(gca,'FontSize',16);
% xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
% ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% % f=gcf;
% % exportgraphics(f,'The_map_zeroAngle_AVG-L.png','Resolution',100)

% %%% plotting include trapping cases 
% figure('color', 'w'); set(gcf, 'Position', [100 100 1200 600]);
% cmap = cmocean('thermal');
% bypass_edge_only = together(:,together(8, :) == 0);
% bypass_tip_only = together(:,together(8, :) == 1);
% scatter(bypass_edge_only(2, :), bypass_edge_only(5, :), 100, bypass_edge_only(3, :), 'Filled','^'); hold on
% scatter(bypass_tip_only(2, :), bypass_tip_only(5, :), 100, bypass_tip_only(3, :), 'Filled'); hold on
% scatter(trapped_together(2, :), trapped_together(5, :), 250, trapped_together(3, :), 'Filled', 'diamond')
% cmap(size(contourL,2)); 
% hcb=colorbar;
% title(hcb,'$\chi_0 (^{\circ})$','FontSize', 16,'Interpreter', 'latex'); % for trapping case plotting.
% grid on
% set(gca,'FontSize',16);
% xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
% ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% % f=gcf;
% % exportgraphics(f,'Chi0-ContourL-InitialY_AVG-L_without-incomplete-trajectory_alldata(till20221005).png','Resolution',100)

%%% the map: contour length: 0.5 < L < 1 & initial position: 0 < y_0 < 1
names(together(2, :) < 0.5) = [];
together(:, together(2, :) < 0.5) = []; % the contour length: 0.5 < L < 1
names(together(2, :) > 1) = [];
together(:, together(2, :) > 1) = [];
names(together(5, :) < 0) = [];
together(:, together(5, :) < 0) = [];  % the initial position: 0 < y_0 < 1
names(together(5, :) > 1) = [];
together(:, together(5, :) > 1) = [];

figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
cmap = cmocean('thermal');
bypass_edge_only = together(:,together(8, :) == 0);
bypass_tip_only = together(:,together(8, :) == 1);
scatter(bypass_edge_only(2, :), bypass_edge_only(5, :), 200, bypass_edge_only(1, :), 'Filled','^'); hold on
scatter(bypass_tip_only(2, :), bypass_tip_only(5, :), 200, bypass_tip_only(1, :), 'Filled'); hold on
cmap(size(contourL,2)); 
hcb=colorbar;
title(hcb,'$Deviation\ (\delta/h_{obs})$','FontSize', 16,'Interpreter', 'latex');
grid on
set(gca,'FontSize',16);
xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,['Deviation-ContourL-InitialY_Angle',num2str(range_chi_up), ...
%     num2str(range_chi_low),'_AVG-L_without-incomplete-trajectory_alldata(till20221005).png'],'Resolution',100)

%%% pick the data point on the last plot to get the case name.
[the_L, the_y0] = ginput(1); % pick up the point you want to show the trajectory.
the_loc = intersect(find(together(2, :)>the_L*0.98 & together(2, :)<the_L*1.02),...
    find(together(5, :)>the_y0*0.98 & together(5, :)<the_y0*1.02));  % Change the control range if there is an error.

names(the_loc)     % show the file

% %%% the speed: contour length: 0.5 < L < 1 & initial position: 0 < y_0 < 1
% figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% cmap = cmocean('thermal');
% scatter(together(2, :), together(5, :), 200, together(9, :), 'Filled','square') % remember to change the 9th row of 'together'
% cmap(size(contourL,2)); 
% hcb=colorbar;
% title(hcb,'$Average\ speed\ (\mu{m}/s)$','FontSize', 16,'Interpreter', 'latex');
% grid on
% set(gca,'FontSize',16);
% xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
% ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% % f=gcf;
% % exportgraphics(f,'Avespeed-ContourL-InitialY_Angle10-10_contourL0.5-1_Y0-1_AVG-L_without-incomplete-trajectory_alldata(till20221005).png','Resolution',100)

%%% 1 - U_bar/U0: contour length: 0.5 < L < 1 & initial position: 0 < y_0 < 1
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
cmap = cmocean('thermal');
scatter(together(2, :), together(5, :), 200, together(9, :), 'Filled','square') % remember to change the 9th row of 'together'
cmap(size(contourL,2)); 
hcb=colorbar;
title(hcb,'$1-\bar{U}/U_0$','FontSize', 16,'Interpreter', 'latex');
grid on
set(gca,'FontSize',16);
xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
f=gcf;
exportgraphics(f,['U_barU0-ContourL-InitialY_Angle',num2str(range_chi_up),num2str(range_chi_low),'' ...
    '_contourL0.5-1_Y0-1_AVG-L_without-incomplete-trajectory_alldata(till20221005).png'],'Resolution',100)


%%% the delta Chi: contour length: 0.5 < L < 1 & initial position: 0 < y_0 < 1
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
cmap = cmocean('thermal');
scatter(together(2, :), together(5, :), 200, together(10, :)', 'Filled','diamond')    % add ' to avoid obfuscating to a triplet.
cmap(size(contourL,2)); 
hcb=colorbar;
title(hcb,'$\chi_f - \chi_0\ (^{\circ})$','FontSize', 16,'Interpreter', 'latex');
grid on
set(gca,'FontSize',16);
xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
f=gcf;
exportgraphics(f,['DeltaChi-ContourL-InitialY_Angle',num2str(range_chi_up),num2str(range_chi_low), ...
    '_contourL0.5-1_Y0-1_AVG-L_without-incomplete-trajectory_alldata(till20221005).png'],'Resolution',100)

%%% the |delta Chi| vs deviation (contour length: 0.5 < L < 1 & initial position: 0 < y_0 < 1)
figure('color', 'w'); set(gcf, 'Position', [100 100 600 400]);
plot(abs(together(10, :)), together(1, :), 'Color','r', 'LineStyle','none', 'Marker','.', 'MarkerSize', 30)
set(gca,'FontSize',16);
xlabel('$\left|\chi_f - \chi_0 \right| (^{\circ})$','FontSize',22,'Interpreter', 'latex');
ylabel('$Deviation\ (\delta/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,['DeltaChi-Deviation_Angle',num2str(range_chi_up),num2str(range_chi_low), ...
%     '_contourL0.5-1_Y0-1_AVG-L_without-incomplete-trajectory_alldata(till20221005).png'],'Resolution',100)

%%% the 1 - U_bar/U0 vs deviation (contour length: 0.5 < L < 1 & initial position: 0 < y_0 < 1)
figure('color', 'w'); set(gcf, 'Position', [100 100 600 400]);
plot(together(9, :), together(1, :), 'Color','r', 'LineStyle','none', 'Marker','.', 'MarkerSize', 30)
set(gca,'FontSize',16);
xlabel('$1-\bar{U}/U_0$','FontSize',22,'Interpreter', 'latex');
ylabel('$Deviation\ (\delta/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
xlim([0 0.6])
% f=gcf;
% exportgraphics(f,['U_barU0-Deviation_Angle',num2str(range_chi_up),num2str(range_chi_low), ...
%     '_contourL0.5-1_Y0-1_AVG-L_without-incomplete-trajectory_alldata(till20221005).png'],'Resolution',100)
