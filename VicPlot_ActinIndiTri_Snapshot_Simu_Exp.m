%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% calculate contact infromation %%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

mag = 0.1; % um/pixel

% colormap
thecolor = colormap('turbo');

% load exp data
[exp_filename, exp_pathname] = uigetfile(['F:\Processing & Results\' ...
    'FSI - Actin &  Individual Obstacle\20230405-Actin-Individual_triangular' ...
    'Pillar_uppoint\results_plus\*.mat'], ...
    'Please choose the file accordingly and carefully!!');
load(fullfile(exp_pathname, exp_filename));

% load image
[tif_filename, tif_pathname] = uigetfile(['F:\Experimental Data (EXTRACTED)' ...
    '\FSI - Actin &  Individual Obstacle\20230405-Actin-Individual_triangular' ...
    'Pillar_uppoint\AfterAveBGR\*.tif'], ...
    'Please choose the tif accordingly and carefully!!');

% load simu data
simu_pathname = uigetdir(['D:\Dropbox\Collaboration - LadHyX\' ...
    'Give_to_Zhibo_nonShared\Data_Give_to_Zhibo_Actin_Indi_Pillar_202306\'], ...
    'Please choose the folder accordingly and carefully!!');
simu_fileinfo = dir(fullfile(simu_pathname, 'output_data\*.vtk'));

% load obstacle data
load(['F:\Processing & Results\FSI - Actin &  Individual Obstacle\' ...
    '20230405-Actin-Individual_triangularPillar_uppoint\' ...
    '20230405_M63_TriUp_1nM_1nL_Expo20ms_1_skel.mat'])
dilated_tri = imdilate(tri/max(max(tri)), strel('disk', 3));
% imshow(dilated_tri, [])

% chosen_exp_case = [3 9 14 18 25 40 137 151 157]; % 20230405_case_51
chosen_exp_case = [2 7 12 17 22 30 41 50 62 72 78]; % 20230405_case_61

n = 1;
for frm_ind = 1:size(Good_case_frm,2)

    xy_ind = Good_case_frm(frm_ind); % index of the 'good' cases
    image_no = xy.frame(xy_ind);
    if ismember(image_no, chosen_exp_case)

        spl = xy.spl{xy_ind}; % the x-y coordinates of the spline (in 'plot' coordinate)
        spl(:, 1) = spl(:, 1); % lzero is very important !!!
        spl(:, 2) = 2048 - spl(:, 2); % convert the x-y coordinates to 'image' coordinate
        mask = zeros(2048, 2048);
        mask(sub2ind([2048 2048], round(spl(:,2)), round(spl(:,1)))) = 1;

        dilatedmask = imdilate(mask, strel('disk', 13));

        if frm_ind == 1
            II_base = imread(fullfile(tif_pathname, tif_filename), xy.frame(xy_ind));
            II_base = double(double(II_base) / (2^16 - 1));
        else
            II = imread(fullfile(tif_pathname, tif_filename), image_no);
            II = double(double(II) / (2^16 - 1));
            II_base = II_base .* (1-dilatedmask) + II .* dilatedmask;
        end

        EXP_center(n, :) = xy.centroid{xy_ind};
        n = n + 1;

    end
end

%% run this part to compare and get 'chosen_simu_case_no'

for ii = 1: length(simu_fileinfo)

    snapshot = readVTK(fullfile(simu_fileinfo(ii).folder, simu_fileinfo(ii).name));
    XY = snapshot.points(:, 1:2)*1e6;
    if ii == 1
        simu_1st_center = mean(XY, 1);
        move_dis = simu_1st_center - EXP_center(1, :)*mag;
    end
    simu_center = mean(XY, 1);
    ShowSimuFiberCoM_x(ii) = simu_center(1) - simu_1st_center(1) + move_dis(1);

end
ShowSimuFiberCoM_x = ShowSimuFiberCoM_x - ShowSimuFiberCoM_x(1);
ShowExpFiberCoM_x = (EXP_center(:,1) - EXP_center(1,1))*mag;

% chosen_simu_case_no = [1 13 23 33 50 109 220 248 262]; % for no_51
chosen_simu_case_no = [1 11 21 31 41 59 84 108 134 155 168]; % for no_61
 
%%

color_max = max(chosen_exp_case(end) * 2, chosen_simu_case_no(end));

figure('color', 'w'); set(gcf, 'Position', [100 300 1000 150]);
n = 1;
EXP_center(:, 2) = 2048 - EXP_center(:, 2);
for ii = chosen_simu_case_no

    color_no_simu = round((chosen_simu_case_no(n) - chosen_simu_case_no(1))/color_max*255) + 1;
    color_no_exp = round((chosen_exp_case(n)*2 - chosen_exp_case(1))/color_max*255) + 1;

    snapshot = readVTK(fullfile(simu_fileinfo(ii).folder, simu_fileinfo(ii).name));
    XY = snapshot.points(:, 1:2)*1e6;
    if ii == 1
        simu_1st_center = mean(XY, 1);
        move_dis = EXP_center(1, :)*mag - simu_1st_center;
    end
    plot(XY(:, 1), XY(:, 2), 'Color', thecolor(color_no_simu, :), 'LineWidth', 2); hold on
    axis equal; axis off
    plot(EXP_center(n, 1)*mag - move_dis(1), EXP_center(n, 2)*mag - move_dis(2), ...
        '.', 'Color', thecolor(color_no_exp, :), 'MarkerSize', 20); hold on
    n = n + 1;

end
% plot obstacle
plot(obs_2d(:,1)*mag - move_dis(1), obs_2d(:,2)*mag - move_dis(2), ':k','LineWidth',1)
% xlim([500 700]); ylim([190 220]);

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['D:\Dropbox\Research\My PhD thesis\' ...
    'Figures\3-rigid_fiber\actin\actin_indiTri_simu_exp_comparison_.eps']);

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['D:\Dropbox\Research\My PhD thesis\Figures\3-rigid_fiber\actin' ...
    '\actin_indiTri_simu_exp_comparison_.pdf']);

%% plot colorbar
figure('color', 'w'); set(gcf, 'Position', [100 100 650 200]);

colormap('turbo');
caxis([0 color_max/100])
c = colorbar;
c.Label.String = 'Time (s)';
c.Label.Interpreter = 'LaTeX';
c.TickLabelInterpreter = 'LaTeX';
c.FontSize = 14;
c.Location = 'west';
axis off

f=gcf;
set(f,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['D:\Dropbox\Research\My PhD thesis\' ...
    'Figures\3-rigid_fiber\actin\actin_indiTri_simu_exp_comparison_verticalColorbar_.eps']);
%%

II_overlap = dilated_tri * (max(max(II_base))/ 10) + II_base;
% figure('color', 'w');  imshow(flipud(II_overlap(850:1150, :)), []); % for no_51
% figure('color', 'w');  imshow(flipud(II_overlap(880:1180, :)), []); % for no_61

% set(gcf,'renderer','Painters');
% print('-depsc2','-tiff','-r100','-vector',['D:\Dropbox\Research\My PhD thesis\' ...
%     'Figures\3-rigid_fiber\actin\actin_indiTri_exp.eps']);

% imwrite(imadjust(flipud(II_overlap(850:1150, :))), ['D:\Dropbox\Research\My PhD thesis\' ...
%     'Figures\3-rigid_fiber\actin\actin_indiTri_exp.tif']); % for no_51
imwrite(imadjust(flipud(II_overlap(880:1180, :))), ['D:\Dropbox\Research\My PhD thesis\' ...
    'Figures\3-rigid_fiber\actin\actin_indiTri_exp_.tif']); % for no_61

