clear; close all; clc;
pathname = uigetdir('D:\Dropbox\PROCESS remotely\Processing & Results', 'Choose a folder');  % input file
Files = dir(fullfile(pathname, 'results', '*.mat'));
bg = imread(fullfile(pathname, 'results', "The Pillar.tif"));
im_BW = imbinarize(bg);
LL = regionprops(im_BW, 'Centroid');
Pillar_CoM  = LL.Centroid;
the_pillar_im = imcomplement(imfill(im_BW, 'holes'));
cmap = cmocean('thermal');
Obj_Mag = 0.63;

im_no = 1:length(Files);
for ii = 1:length(Files)

    load([Files(ii).folder, filesep, Files(ii).name]);

    %%% Plot trajectories
%     figure('color', 'w'); set(gcf, 'Position', [300 300 600 100]);
%     imshow(the_pillar_im(801:1200, :)); hold on
% 
%     for k =1:5:min(128, xy.nframe) % 128: because of the colormap
%         plot(xy.spl{k}(:,1),2048-xy.spl{k}(:,2)-800,'Color',cmap(k*2,:),'LineWidth',1.5)
%         hold on
%     end
%     axis equal
%     hold off
% 
%     f=gcf;
%     exportgraphics(f,['E:\Dropbox\PROCESS remotely\Processing\20220624-SU8_Fibers-' ...
%         'Individual_triangularPillar_uppoint\results\',Files(ii).name,'.png'],'Resolution',1000)
%     close

    All_data.filename{ii} = Files(ii).name;    

    sorted_lengths = sort(xy.arclen_spl(Good_case_frm));
    All_data.contour_length(ii) = max(sorted_lengths) * Obj_Mag;  % Select the lengest filaments as the contour length (UNIT: um).

    centroidxy = reshape(cell2mat(xy.centroid),2,numel(xy.centroid));
    centroidxy = centroidxy(:, Good_case_frm);
    All_data.delta_y(ii) = abs(centroidxy(2, 1) - centroidxy(2,end)) * Obj_Mag;
    All_data.initial_y(ii) = ((2048-centroidxy(2, 1)) - Pillar_CoM(2)) * Obj_Mag + 25; % Definition of y0 is the distance between CoM and the edge, so +25.
    All_data.CoM{ii} = centroidxy;
    %     hold on; plot(centroidxy(1,:), 2048-centroidxy(2,:)-800, "Color",'k', 'LineStyle','--');

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

% save(fullfile(pathname, 'results', 'infos.mat'), 'All_data');

InfoMatrix = [im_no; All_data.contour_length; All_data.delta_y; All_data.initial_y; Chi];
[temp, order] = sort(InfoMatrix(3,:));
sorted_InfoMatrix = InfoMatrix(:,order);

figure('color', 'w');
plot(sorted_InfoMatrix(2,:), sorted_InfoMatrix(3,:), "LineStyle", "none", "Marker",".", 'MarkerSize', 30)
set(gca,'FontSize',14);
xlabel('$Contour\ length\ L\ ({\mu}m)$','FontSize', 20,'Interpreter', 'latex');
ylabel('$Deviation\ |y|\ ({\mu}m)$','FontSize', 20,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,['D:\Dropbox\PROCESS remotely\Processing\20220624-SU8_Fibers-' ...
%     'Individual_triangularPillar_uppoint\results\L_deltaY.png'],'Resolution',100);
% close

figure('color', 'w');
plot(sorted_InfoMatrix(4,:), sorted_InfoMatrix(3,:), "LineStyle", "none", "Marker",".", 'MarkerSize', 30)
xlabel('$Initial\ position\ y_0\ ({\mu}m)$','FontSize', 20,'Interpreter', 'latex');
ylabel('$Deviation\ |y|\ ({\mu}m)$','FontSize', 20,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,['D:\Dropbox\PROCESS remotely\Processing\20220624-SU8_Fibers-' ...
%     'Individual_triangularPillar_uppoint\results\Y0_deltaY.png'],'Resolution',100)
% close

figure('color', 'w');
plot(sorted_InfoMatrix(5,:), sorted_InfoMatrix(3,:), "LineStyle", "none", "Marker",".", 'MarkerSize', 30)
xlabel('$Initial\ angle\ \chi_0\ (^{\circ})$','FontSize', 20,'Interpreter', 'latex');
ylabel('$Deviation\ |y|\ ({\mu}m)$','FontSize', 20,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,['D:\Dropbox\PROCESS remotely\Processing\20220624-SU8_Fibers-' ...
%     'Individual_triangularPillar_uppoint\results\Chi0_deltaY.png'],'Resolution',100)
% close


sorted_InfoMatrix(:, (sorted_InfoMatrix(5,:) > 10)) = []; sorted_InfoMatrix(:, (sorted_InfoMatrix(5,:) < -10)) = [];
sorted_InfoMatrix(:, (sorted_InfoMatrix(4,:) > 12.5)) = []; sorted_InfoMatrix(:, (sorted_InfoMatrix(4,:) < -6.25)) = [];
figure('color', 'w');
plot(sorted_InfoMatrix(2,:), sorted_InfoMatrix(3,:), "LineStyle", "none", "Marker",".", 'MarkerSize', 30, 'Color','r')
set(gca,'FontSize',16); set(gcf, 'Position', [100 100 800 600]);
title('$-10^{\circ} < \chi_0 < 10^{\circ}\ and\ 0.25 < y_0 / h_{obs} < 0.5$','FontSize', 20,'Interpreter', 'latex')
xlabel('$Contour\ length\ L\ ({\mu}m)$','FontSize', 20,'Interpreter', 'latex');
ylabel('$Deviation\ |y|\ ({\mu}m)$','FontSize', 20,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,['G:\PhD, PMMH, ESPCI\Processing\20220624-SU8_Fibers-Individual_' ...
%     'triangularPillar_uppoint\results\L_deltaY_selected_Y0basedontheedge.png'],'Resolution',100)
% close


sorted_InfoMatrix(:, (sorted_InfoMatrix(5,:) > 10)) = []; sorted_InfoMatrix(:, (sorted_InfoMatrix(5,:) < -10)) = [];
sorted_InfoMatrix(:, (sorted_InfoMatrix(2,:) < 80)) = []; sorted_InfoMatrix(:, (sorted_InfoMatrix(2,:) > 120)) = [];
figure('color', 'w');
% plot the normalized y0 (from the edge of the triangle pillar) based on
% the h_obs (= 75um)
plot((sorted_InfoMatrix(4,:)+25)/75, sorted_InfoMatrix(3,:), "LineStyle", "none", "Marker",".", 'MarkerSize', 30, 'Color','r')
set(gca,'FontSize',16); set(gcf, 'Position', [100 100 800 600]);
title('$-10^{\circ} < \chi_0 < 10^{\circ}\ and\ 80{\mu}m < L < 120{\mu}m$','FontSize', 20,'Interpreter', 'latex')
xlabel('$Initial\ position\ y_0 / h_{obs}$','FontSize', 20,'Interpreter', 'latex');
ylabel('$Deviation\ |y|\ ({\mu}m)$','FontSize', 20,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,['G:\PhD, PMMH, ESPCI\Processing\20220624-SU8_Fibers-Individual_' ...
%     'triangularPillar_uppoint\results\Y0_basedontheedge_deltaY_selected.png'],'Resolution',100)
% close

%% Drawing
clear; close all; clc;

h_obs = 75; l_obs = 86.6;

load('20220624-SU8_Fibers-Individual_triangularPillar_uppoint\results\infos.mat');
All_data1 = All_data;
load('20220913-SU8_Fibers-Individual_triangularPillar_uppoint\results\infos.mat');
All_data2 = All_data;
All_data = [All_data1, All_data2];

contourL = [All_data(1).contour_length, All_data(2).contour_length];
norm_contourL = contourL / l_obs;
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
edges = 20:20:180;
histogram(contourL,edges);
set(gca,'FontSize',16);
xlabel('$Contour\ length\ L\ ({\mu}m)$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Number\ of\ cases$','FontSize', 22,'Interpreter', 'latex');
xlim([20 180]);
f=gcf;
exportgraphics(f,'Statistics_contourL.png','Resolution',100)

Chi = [All_data(1).Chi, All_data(2).Chi];
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
edges = -90:20:90;
histogram(Chi,edges);
set(gca,'FontSize',16);
xlabel('$Initial\ angle\ \chi_0\ (^{\circ})$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Number\ of\ cases$','FontSize', 22,'Interpreter', 'latex');
xlim([-90 90]);
f=gcf;
exportgraphics(f,'Statistics_Chi0.png','Resolution',100)

initial_y = [All_data(1).initial_y, All_data(2).initial_y];
norm_initial_y = initial_y / h_obs;
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
edges = -0.2:0.2:1.2;
histogram(norm_initial_y,edges);
set(gca,'FontSize',16);
xlabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Number\ of\ cases$','FontSize', 22,'Interpreter', 'latex');
xlim([-0.2 1.2]);
f=gcf;
exportgraphics(f,'Statistics_y0.png','Resolution',100)

delta_y = [All_data(1).delta_y, All_data(2).delta_y];
norm_delta_y = delta_y / h_obs;

% the map: without any 'filter'.
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
cmap = cmocean('thermal');
scatter(norm_contourL, norm_initial_y, 50, norm_delta_y, 'Filled')
cmap(size(contourL,2)); 
hcb=colorbar;
title(hcb,'$Deviation\ (\delta/h_{obs})$','FontSize', 16,'Interpreter', 'latex');
grid on
set(gca,'FontSize',16);
xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
f=gcf;
exportgraphics(f,'The_map_full.png','Resolution',100)


% the map: -10 < Chi < 10.
together = [norm_contourL; norm_initial_y; Chi; norm_delta_y];
together(:, together(3, :) < -10) = [];
together(:, together(3, :) > 10) = [];
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
cmap = cmocean('thermal');
scatter(together(1, :), together(2, :), 50, together(4, :), 'Filled')
cmap(size(contourL,2)); 
hcb=colorbar;
title(hcb,'$Deviation\ (\delta/h_{obs})$','FontSize', 16,'Interpreter', 'latex');
grid on
set(gca,'FontSize',16);
xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
f=gcf;
exportgraphics(f,'The_map_zeroAngle.png','Resolution',100)

% the map -10 < Chi < 10.
together = [norm_contourL; norm_initial_y; Chi; norm_delta_y];
together(:, together(3, :) < -10) = [];
together(:, together(3, :) > 10) = [];
together(:, together(1, :) < 0.5) = [];
together(:, together(1, :) > 1) = [];
together(:, together(2, :) < 0) = [];
together(:, together(2, :) > 1) = [];
figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
cmap = cmocean('thermal');
scatter(together(1, :), together(2, :), 50, together(4, :), 'Filled')
cmap(size(contourL,2)); 
hcb=colorbar;
title(hcb,'$Deviation\ (\delta/h_{obs})$','FontSize', 16,'Interpreter', 'latex');
grid on
set(gca,'FontSize',16);
xlabel('$Contour\ length\ (L/l_{obs})$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Initial\ position\ (y_0/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
f=gcf;
exportgraphics(f,'The_map_zeroAngle_Simu.png','Resolution',100)


