%%%% Poincaré plot or recurrence plot (RP).

clear; close all; clc;

xlsfile = readcell('ForActinPostprocessing.xlsx','Sheet','Sheet1','NumHeaderLines',1);
% This is the file that contains all the information about the later processing (in sheet 1).

NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.
PAsPath = xlsfile(:, 3);  % Path of the pillar array information.

n = 1;
for no_Group = [7 8 14 15 16 17 18]
    % Square-based array 0°, 10°, and 20°
    % No.13 needed to be recalculated (20230125)

    if no_Group == 7 || no_Group == 8
        Array_angle = 0;
    elseif no_Group == 14 || no_Group == 15
        Array_angle = 10;
    elseif no_Group == 16 || no_Group == 17 || no_Group == 18
        Array_angle = 20;
    end
    RotMatrix_init = rotz(-Array_angle); RotMatrix_init = RotMatrix_init(1:2, 1:2);
    % to rotate the pillar array

    the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
    filelist = dir(fullfile(storePath{no_Group},'*.mat'));
    % list of the .mat files which contain the reconstruction information
    % (came from 'Filaments detection' code) in one group.

    for no_Case = 1:length(filelist)

        filename = filelist(no_Case).name;
        if contains(filename, 'PAsInfoAdded_')
            load([storePath{no_Group}, filesep , filename])
            centers(:, 2) = 2048 - centers(:, 2);  % flip to image coordinate
            centers_new = (RotMatrix_init * centers')'; % rotate based on the design
            % viscircles(centers_new, radii,'LineStyle','--', 'LineWidth', 0.5, 'Color', 'k'); axis equal

            % correct the rotation degree
            sorted_centers_new = sortrows(centers_new,2);
            sorty = sort(centers_new(:, 2));
            sorty_ind = find(diff(sorty) > 100);
            corrected_angle = zeros(length(sorty_ind)-1, 1);
            for ii = 1: length(sorty_ind)-1
                to_be_fitted = sorted_centers_new(sorty_ind(ii)+1:sorty_ind(ii+1), :);
                fit_linear = fit(to_be_fitted(:, 1), to_be_fitted(:, 2), 'poly1');
                k = fit_linear.p1;
                corrected_angle(ii) = atand(k);
            end
            corrected_angle = mean(corrected_angle);
            RotMatrix_correct = rotz(-Array_angle-corrected_angle);
            RotMatrix_correct = RotMatrix_correct(1:2, 1:2);
            centers_corrected = (RotMatrix_correct * centers')';
            % viscircles(centers_corrected, radii,'LineStyle','--', 'LineWidth', 0.5, 'Color', 'r'); axis equal

            % get the approximate (after-gridlization) pillar center positions
            sortxx = sort(centers_corrected(:, 1)); sortyy = sort(centers_corrected(:, 2));
            sortxx_ind = find(diff(sortxx) > 100); sortyy_ind = find(diff(sortyy) > 100);
            PAs_X = zeros(length(sortxx_ind)+1, 1); PAs_Y = zeros(length(sortyy_ind)+1, 1);
            for jj = 1: length(sortxx_ind)-1
                PAs_X(jj+1) = mean(sortxx(sortxx_ind(jj)+1:sortxx_ind(jj+1)));
            end
            for jj = 1: length(sortyy_ind)-1
                PAs_Y(jj+1) = mean(sortyy(sortyy_ind(jj)+1:sortyy_ind(jj+1)));
            end
            PAs_X(1) = PAs_X(2) - mean(diff(PAs_X(2:end-1)));
            PAs_X(end) = PAs_X(end-1) + mean(diff(PAs_X(2:end-1)));
            PAs_Y(1) = PAs_Y(2) - mean(diff(PAs_Y(2:end-1)));
            PAs_Y(end) = PAs_Y(end-1) + mean(diff(PAs_Y(2:end-1)));
            % [XX, YY] = meshgrid(PAs_X, PAs_Y); plot(XX(:), YY(:))

            % average gap along y-direction
            ave_y_gap = mean(diff(PAs_Y));

            % get the filament length
            spl_Ls = xy.arclen_spl(Good_case_frm);
            L_0 = Vic_Get_ave_cutExtreme(spl_Ls, 0.2);

            % get the filament CoMs at different time
            fiber_center = reshape(cell2mat(xy.centroid), 2, length(cell2mat(xy.centroid))/2)';
            fiber_center = fiber_center(Good_case_frm, :);
            fiber_center = fiber_center + [lzero, lzero];
            fiber_center = (RotMatrix_correct * fiber_center')';
            % plot(fiber_center(:, 1), fiber_center(:, 2)); hold on
            % viscircles(centers_corrected, radii,'LineStyle','--', 'LineWidth', 0.5, 'Color', 'r'); axis equal

            % divide the trajectory into pieces according to their y-position
            fiber_Y_indicator = fiber_center(:, 2);
            for kk = 1: length(PAs_Y)-1
                fiber_Y_indicator(fiber_Y_indicator > PAs_Y(kk) & fiber_Y_indicator < PAs_Y(kk+1)) = kk;
            end

            % the 'entering lattice' positions
            ind_ToBeMoved = min(abs(repmat(fiber_center(:, 1), 1, length(PAs_X)) - PAs_X'), [], 2) > 60;
            fiber_center(ind_ToBeMoved, :) = [];  % only keep the cases that close to the lattice verticle edge
            fiber_Y_indicator(ind_ToBeMoved) = [];  % Y indicator as well
            fiber_X = fiber_center(:, 1);
            For_Poincare = nan(length(PAs_X), 3);
            for kk = 1: length(PAs_X)
                to_be_fitted2 = fiber_center(abs(fiber_X-PAs_X(kk))<=60, :);
                if ~isempty(to_be_fitted2) && numel(to_be_fitted2) > 2
                    fit_linear2 = fit(to_be_fitted2(:, 1), to_be_fitted2(:, 2), 'poly1');
                    For_Poincare(kk, 1) = mod(fit_linear2(PAs_X(kk)), ave_y_gap) / ave_y_gap;
                    For_Poincare(kk, 2) = kk;
                    For_Poincare(kk, 3) = fiber_Y_indicator(kk);
                end
            end
            For_Poincare(:, 3) = For_Poincare(:, 3) - min(For_Poincare(:, 3)) + 1;

            Info.L(n) = L_0;
            Info.map{n} = For_Poincare;
            Info.date(n) = the_exp_date;
            n = n + 1;

            % f=gcf;
            % exportgraphics(f,['D:\Dropbox\GitHub\tmp', filesep, num2str(the_exp_date), filename(25:end-11), '.png'],'Resolution',100)

        end
    end
end

save('Poincare_Map_data.mat', 'Info')

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