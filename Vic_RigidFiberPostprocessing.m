clear; close all; clc;
pathname = uigetdir('F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle', 'Choose a folder');  % input file
Files = dir(fullfile(pathname, 'results with timestamps', '*.mat'));
bg = imread(fullfile(pathname, 'results', "The Pillar.tif"));
im_BW = imbinarize(bg);
LL = regionprops(im_BW, 'Centroid');
Pillar_CoM  = LL.Centroid;
the_pillar_im = imcomplement(imfill(im_BW, 'holes'));
% cmap = cmocean('thermal');
Obj_Mag = 0.63;

im_no = 1:length(Files);
for ii = 1:length(Files)

    load([Files(ii).folder, filesep, Files(ii).name]);

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

    All_data.filename{ii} = Files(ii).name;    

    sorted_lengths = sort(xy.arclen_spl(Good_case_frm));
    All_data.contour_length(ii) = mean(sorted_lengths) * Obj_Mag;  % Average to get the contour length (UNIT: um).

    centroidxy = reshape(cell2mat(xy.centroid),2,numel(xy.centroid));
    centroidxy = centroidxy(:, Good_case_frm);
    All_data.delta_y(ii) = abs(centroidxy(2, 1) - centroidxy(2,end)) * Obj_Mag;
    All_data.initial_y(ii) = ((2048-centroidxy(2, 1)) - Pillar_CoM(2)) * Obj_Mag + 25; % Definition of y0 is the distance between CoM and the edge, so +25.
    All_data.initial_x(ii) = (centroidxy(1, 1) - Pillar_CoM(1)) * Obj_Mag;
    All_data.final_y(ii) = ((2048-centroidxy(2, end)) - Pillar_CoM(2)) * Obj_Mag + 25; % Definition of y0 is the distance between CoM and the edge, so +25.
    All_data.final_x(ii) = (centroidxy(1, end) - Pillar_CoM(1)) * Obj_Mag;
    All_data.CoM{ii} = centroidxy;
    %     hold on; plot(centroidxy(1,:), 2048-centroidxy(2,:)-800, "Color",'k', 'LineStyle','--');

    % check if the fiber bypass the obtacle tip
    if 2048 - mean(centroidxy(2, and(Pillar_CoM(1)-100 < centroidxy(1, :), centroidxy(1, :) < Pillar_CoM(1)+100))) > Pillar_CoM(2) % if it goes below (in image system).
        bypass_tip = true;
    else
        bypass_tip = false;
    end
    All_data.bypass_tip(ii) = bypass_tip;

    crd = xy.crd{1,1};
    centroidxy_k = xy.centroid{1,1};
    Gyr = 1/size(crd,1) * [sum((crd(:, 1)-centroidxy_k(1)).^2),  sum((crd(:, 1)-centroidxy_k(1)) .* (crd(:, 2)-centroidxy_k(2)));
        sum((crd(:, 2)-centroidxy_k(2)) .* (crd(:, 1)-centroidxy_k(1))), sum((crd(:, 2)-centroidxy_k(2)).^2)];

    [eigenV,eigenD] = eig(Gyr);
    [d,ind] = sort(diag(eigenD));
    Ds = eigenD(ind,ind);
    Vs = eigenV(:,ind);
    All_data.Chi(ii)  = atan(Vs(2,2)/Vs(1,2))/pi*180;

end

% save([pathname, '_infos_contourL-fromAVG.mat'], 'All_data');



%% Ploting
clear; close all; clc;

h_obs = 75; l_obs = 86.6;

load('20220624-SU8_Fibers-Individual_triangularPillar_uppoint_infos_contourL-fromAVG.mat');
All_data1 = All_data;
load('20220913-SU8_Fibers-Individual_triangularPillar_uppoint_infos_contourL-fromAVG.mat');
All_data2 = All_data;
All_data = [All_data1, All_data2];

names = [All_data(1).filename, All_data(2).filename];

delta_y = [All_data(1).delta_y, All_data(2).delta_y];
norm_delta_y = delta_y / h_obs; % normalized delta
contourL = [All_data(1).contour_length, All_data(2).contour_length];
norm_contourL = contourL / l_obs; % normalized L
Chi = [All_data(1).Chi, All_data(2).Chi]; 

initial_x = [All_data(1).initial_x, All_data(2).initial_x];
norm_initial_x = initial_x / h_obs; % normalized x_0
initial_y = [All_data(1).initial_y, All_data(2).initial_y];
norm_initial_y = initial_y / h_obs; % normalized y_0
final_x = [All_data(1).final_x, All_data(2).final_x];
norm_final_x = final_x / h_obs; % normalized x_f
final_y = [All_data(1).final_y, All_data(2).final_y];
norm_final_y = final_y / h_obs; % normalized y_f
bypass_tip = [All_data(1).bypass_tip, All_data(2).bypass_tip];

%%% statistics
% figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% edges = 20:20:180;
% histogram(contourL,edges);
% set(gca,'FontSize',16);
% xlabel('$Contour\ length\ L\ ({\mu}m)$','FontSize', 22,'Interpreter', 'latex');
% ylabel('$Number\ of\ cases$','FontSize', 22,'Interpreter', 'latex');
% xlim([20 180]);
% % f=gcf;
% % exportgraphics(f,'Statistics_contourL.png','Resolution',100)

%%% statistics
% figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% edges = -90:20:90;
% histogram(Chi,edges);
% set(gca,'FontSize',16);
% xlabel('$Initial\ angle\ \chi_0\ (^{\circ})$','FontSize', 22,'Interpreter', 'latex');
% ylabel('$Number\ of\ cases$','FontSize', 22,'Interpreter', 'latex');
% xlim([-90 90]);
% % f=gcf;
% % exportgraphics(f,'Statistics_Chi0.png','Resolution',100)

%%% statistics
% figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
% edges = -0.2:0.2:1.2;
% histogram(norm_initial_y,edges);
% set(gca,'FontSize',16);
% xlabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% ylabel('$Number\ of\ cases$','FontSize', 22,'Interpreter', 'latex');
% xlim([-0.2 1.2]);
% % f=gcf;
% % exportgraphics(f,'Statistics_y0.png','Resolution',100)

together = [norm_delta_y; norm_contourL; Chi; norm_initial_x; norm_initial_y; norm_final_x; norm_final_y; bypass_tip];
names(together(4, :) > -7) = [];
together(:, together(4, :) > -7) = []; % remove the cases too close to the obstacle on the upstream side.
names(together(6, :) < 7) = [];
together(:, together(6, :) < 7) = []; % remove the cases too close to the obstacle on the downstream side.

% the map: without any 'filter'.
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

names(together(3, :) < -10) = [];
together(:, together(3, :) < -10) = [];  % the initial angle: -10 < Chi < 10.
names(together(3, :) > 10) = [];
together(:, together(3, :) > 10) = [];
% the map -10 < Chi < 10.
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

% the map: contour length: 0.5 < L < 1 & initial position: 0 < y_0 < 1
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
scatter(bypass_tip_only(2, :), bypass_tip_only(5, :), 200, bypass_tip_only(1, :), 'Filled')
cmap(size(contourL,2)); 
hcb=colorbar;
title(hcb,'$Deviation\ (\delta/h_{obs})$','FontSize', 16,'Interpreter', 'latex');
grid on
set(gca,'FontSize',16);
xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,'The_map_zeroAngle_Simu_AVG-L_separated.png','Resolution',100)

% pick the data point on the last plot to get the case name.
[the_L, the_y0] = ginput(1); % pick up the point you want to show the trajectory.
the_loc = intersect(find(together(2, :)>the_L*0.98 & together(2, :)<the_L*1.02),...
    find(together(5, :)>the_y0*0.98 & together(5, :)<the_y0*1.02));  % Change the control range if there is an error.

names(the_loc)     % show the file

