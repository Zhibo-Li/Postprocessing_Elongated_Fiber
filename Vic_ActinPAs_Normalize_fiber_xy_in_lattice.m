%%%% Normalize the fiber x-y coordinates in the lattice.
%
% Run this code after Vic_ActinPAs_Normalize_CoM_xy.m
%
% data from PlusInfo_trajectory_..._batch1.mat
% saving name format: PlusInfo_trajectory_..._batch1.mat

clear; close all; clc;

B = 6.9e-26;  % Bending rigidity

xlsfile = readcell('ForActinPostprocessing.xlsx','Sheet','Sheet1','NumHeaderLines',1);
% This is the file that contains all the information about the later processing (in sheet 1).

NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.
Array_angles = xlsfile(:, 14);  % The flow angles.

for no_Group = [7 8 13:28]

    the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
    thefiles = dir(fullfile(storePath{no_Group},'*.mat'));

    Array_angle = Array_angles{no_Group};
    RotMatrix_init = rotz(-Array_angle); RotMatrix_init = RotMatrix_init(1:2, 1:2);

    for file_ind = 1:length(thefiles)

        filename = thefiles(file_ind).name;

        if contains(filename, 'PAsInfoAdded_')

            pathname = thefiles(1).folder;
            save_pathname = strrep(pathname,'results','results_plus');
            filename = thefiles(file_ind).name
            load([save_pathname, filesep, 'PlusInfo_', filename(14: end-4), '.mat'])

            centers_flip = centers;
            centers_flip(:, 2) = 2048 - centers_flip(:, 2);  % flip to image coordinate
            centers_new = (RotMatrix_init * centers_flip')'; % rotate based on the design
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
            centers_corrected = (RotMatrix_correct * centers_flip')';
            % viscircles(centers_corrected, radii,'LineStyle','--','LineWidth', 0.5, 'Color', 'r'); axis equal; hold on

            % get the approximate (after-gridlization) pillar center positions
            sortxx = sort(centers_corrected(:, 1)); sortyy = sort(centers_corrected(:, 2));
            sortxx_ind = find(diff(sortxx) > 100); sortyy_ind = find(diff(sortyy) > 100);
            PAs_X = zeros(length(sortxx_ind)+3, 1); PAs_Y = zeros(length(sortyy_ind)+3, 1);
            for jj = 1: length(sortxx_ind)-1
                PAs_X(jj+2) = mean(sortxx(sortxx_ind(jj)+1:sortxx_ind(jj+1)));
            end
            for jj = 1: length(sortyy_ind)-1
                PAs_Y(jj+2) = mean(sortyy(sortyy_ind(jj)+1:sortyy_ind(jj+1)));
            end
            % make a grid of pillar CoMs (enlarge the PAs area)
            Ctr2Ctr_x = mean(diff(PAs_X(3:end-2))); Ctr2Ctr_y = mean(diff(PAs_Y(3:end-2)));
            PAs_X(2) = PAs_X(3) - Ctr2Ctr_x; PAs_X(1) = PAs_X(2) - Ctr2Ctr_x;
            PAs_X(end-1) = PAs_X(end-2) + Ctr2Ctr_x; PAs_X(end) = PAs_X(end-1) + Ctr2Ctr_x;
            PAs_Y(2) = PAs_Y(3) - Ctr2Ctr_y; PAs_Y(1) = PAs_Y(2) - Ctr2Ctr_y;
            PAs_Y(end-1) = PAs_Y(end-2) + Ctr2Ctr_y; PAs_Y(end) = PAs_Y(end-1) + Ctr2Ctr_y;

            [PAs_X_mesh,PAs_Y_mesh] = meshgrid(PAs_X,PAs_Y);
            %             PAs_xy = complex(PAs_X_mesh(:), PAs_Y_mesh(:));
            PAs_xy = [PAs_X_mesh(:), PAs_Y_mesh(:)];

            for frm_ind = 1:size(Good_case_frm,2)

                xy_ind = Good_case_frm(frm_ind); % index of the 'good' cases
                spl = xy.spl{1,xy_ind};
                Rotated_spl = RotMatrix_correct * spl';
                %                 Rotated_spl = complex(Rotated_spl(1, :), Rotated_spl(2, :));
                Rotated_spl = [Rotated_spl(1, :); Rotated_spl(2, :)]';

                D = pdist2(PAs_xy, Rotated_spl); % distances between pillar and x-y coordinates of actin.

                PAs_x_repeat = repmat(PAs_X_mesh(:), 1, size(Rotated_spl, 1));
                PAs_y_repeat = repmat(PAs_Y_mesh(:), 1, size(Rotated_spl, 1));

                spl_x_repeat = repmat(Rotated_spl(:, 1)', size(PAs_xy, 1), 1);
                spl_y_repeat = repmat(Rotated_spl(:, 2)', size(PAs_xy, 1), 1);

                dist_x = spl_x_repeat - PAs_x_repeat;
                dist_y = spl_y_repeat - PAs_y_repeat;

                real_sign = dist_x > 0;
                imag_sign = dist_y > 0;

                the_dist_sign = (and(real_sign, imag_sign) - 0.5) * 2;
                dist = the_dist_sign .* D; % distances with sign (only lower-left pillar is positive)

                dist(dist < 0) = nan;
                [ref_x_ind, ref_y_ind] = find(dist == min(dist)); % find the closest lower-left pillar index.

                fiber_xy_in_lattice = [diag(dist_x(ref_x_ind, ref_y_ind))/Ctr2Ctr_x, ...
                    diag(dist_y(ref_x_ind, ref_y_ind))/Ctr2Ctr_y];
                xy.fiber_xy_in_lattice{1,xy_ind} = fiber_xy_in_lattice; % save in xy struct

            end

            save([save_pathname, filesep, 'PlusInfo_', filename(14: end-4), '.mat'],'centers', ...
                'circleMask', 'CoM_x', 'Energy', 'Good_case_frm', 'L_ee_norm', ...
                'L_ee_norm_belowOne', 'lzero','metric','radii','xy', 'Chi', 'aniso', ...
                'Rotated_fiber_CoM_xy_normalized', 'RotMatrix_correct', 'Ctr2Ctr_x', 'Ctr2Ctr_y');

            clearvars CoM_x Energy xy centers radii L_ee_norm_belowOne L_ee_norm Chi aniso Rotated_fiber_CoM_xy_normalized
        end
    end
end
