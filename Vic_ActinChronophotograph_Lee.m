%%%% Plot the chronophotograoh and L_ee of actin filaments.
%%%% data from Vic_ActinAddInformationSave.m
%%%% data name format: PAsInfoAdded_trajectory_M63_..._batch1.mat

clear; close all; clc;

xlsfile = readcell('ForActinPostprocessing.xlsx','Sheet','Sheet1','NumHeaderLines',1);
% This is the file that contains all the information about the later processing (in sheet 1).

NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.

for no_Group = [7 8 13 14 15 16 17 18]

    the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
    thefiles = dir(fullfile(storePath{no_Group},'*.mat'));

    for file_ind = 1:length(thefiles)

        filename = thefiles(file_ind).name;

        if contains(filename, 'PAsInfoAdded_')

            load(fullfile(thefiles(1).folder, thefiles(file_ind).name));

            pathname = thefiles(1).folder;
            filename = thefiles(file_ind).name

            figure('color', 'w'); set(gcf, 'Position', [100 300 1600 300]);
            cmap = colormap('jet');

            spl_Ls = xy.arclen_spl(Good_case_frm);
            L_0 = Vic_Get_ave_cutExtreme(spl_Ls, 0.2); % get the filament length

            CoM_x = zeros(size(Good_case_frm,2), 1); L_ee = zeros(size(Good_case_frm,2), 1);
            for frm_ind = 1:size(Good_case_frm,2)

                xy_ind = Good_case_frm(frm_ind);% index of the 'good' cases

                %%%%%%%%%%%%%% plot fiber snapshots %%%%%%%%%%%%%%%%%%%%%%
                if frm_ind == 1
                    plot(xy.spl{xy_ind}(:,1)+lzero(frm_ind),xy.spl{xy_ind}(:,2)+lzero(frm_ind), 'LineWidth', 2);
                    addaxislabel(1,'y (pixel)');
                elseif mod(frm_ind, 3) == 0
                    addaxisplot(xy.spl{xy_ind}(:,1)+lzero(frm_ind),xy.spl{xy_ind}(:,2)+lzero(frm_ind), 1, 'color', cmap(mod(frm_ind*32, 255)+1, :), 'LineWidth', 2);
                end
                hold on

                %%%%%%%%%%%%%% get L_ee and x-position %%%%%%%%%%%%%%%%%%%
                CoM_x(frm_ind) = xy.centroid{1, xy_ind}(1); % x-position of fiber center-of-mass
                L_ee(frm_ind) = sqrt((xy.spl{xy_ind}(1,1)-xy.spl{xy_ind}(end,1))^2 ...
                    + (xy.spl{xy_ind}(1,2)-xy.spl{xy_ind}(end,2))^2) / L_0;
                if L_ee(frm_ind) > 1
                    L_ee(frm_ind) = 1;
                end

            end
            xlim([0 2050])
            xlabel('x (pixel)')
            axis equal; hold on

            addaxis(CoM_x+lzero(frm_ind), L_ee, [0 1],'*r', 'LineStyle','none', 'MarkerSize', 7);
            addaxislabel(2, 'L_e_e / L_0');
            centers(:, 2) = 2048 - centers(:, 2);
            viscircles(centers, radii,'LineStyle','--', 'LineWidth', 0.5, 'Color', 'k'); hold on

            f=gcf;
            exportgraphics(f,[pathname, filesep, filename(14: end-4), '_Lee_PossibleGreaterThanOne.png'],'Resolution',100)

            close all
            clearvars CoM_x Energy xy centers radii
        end
    end
end





function L_0 = Vic_Get_ave_cutExtreme(spl_Ls, threshold)
% the function is used to calculte the contour lenght of the filament based
% on the average without the extrame values.
%
% threshold -- threshold to truncate the extreme values

L_0_coarse = mean(spl_Ls, 'omitnan'); % roughly estimated filament length
spl_Ls(spl_Ls < L_0_coarse*(1-threshold)) = [];
spl_Ls(spl_Ls > L_0_coarse*(1+threshold)) = []; % remove the extreme value
if isempty(spl_Ls)
    L_0 = L_0_coarse;
    return
end
L_0 = mean(spl_Ls, 'omitnan');

while L_0_coarse ~= L_0 && ~isnan(L_0)
    L_0_coarse = L_0;
    spl_Ls(spl_Ls < L_0_coarse*(1-threshold)) = [];
    spl_Ls(spl_Ls > L_0_coarse*(1+threshold)) = []; % remove the extreme value
    if isempty(spl_Ls)
        L_0 = L_0_coarse;
        return
    end
    L_0 = mean(spl_Ls, 'omitnan');
end
end