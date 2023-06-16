%%%% Normalize the positions based on the pillar array.
%%%% 
%
% data from PlusInfo_trajectory_..._batch1.mat
% saving name format: PlusInfo_trajectory_..._batch1.mat
%

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
            Ctr2Ctr_x = diff(PAs_X(3:end-2)); Ctr2Ctr_y = diff(PAs_Y(3:end-2));
            PAs_X(2) = PAs_X(3) - mean(Ctr2Ctr_x); PAs_X(1) = PAs_X(2) - mean(Ctr2Ctr_x);
            PAs_X(end-1) = PAs_X(end-2) + mean(Ctr2Ctr_x); PAs_X(end) = PAs_X(end-1) + mean(Ctr2Ctr_x);
            PAs_Y(2) = PAs_Y(3) - mean(Ctr2Ctr_y); PAs_Y(1) = PAs_Y(2) - mean(Ctr2Ctr_y);
            PAs_Y(end-1) = PAs_Y(end-2) + mean(Ctr2Ctr_y); PAs_Y(end) = PAs_Y(end-1) + mean(Ctr2Ctr_y);

            fiber_CoM_xy = reshape(cell2mat(xy.centroid(1,Good_case_frm)),2,[]); 
            Rotated_fiber_CoM_xy = RotMatrix_correct * fiber_CoM_xy; % fiber CoM after rotation

            [PAs_X_mesh,PAs_Y_mesh] = meshgrid(PAs_X,PAs_Y);
            PAs_xy = complex(PAs_X_mesh(:), PAs_Y_mesh(:));
            Rotated_fiber_CoM_xy_1st = complex(Rotated_fiber_CoM_xy(1, 1),Rotated_fiber_CoM_xy(2, 1));

            % find the lower-left, closest pillar CoM of the first fiber CoM
            real_sign = real(Rotated_fiber_CoM_xy_1st - PAs_xy) > 0;
            imag_sign = imag(Rotated_fiber_CoM_xy_1st - PAs_xy) > 0;
            the_dist_sign = (and(real_sign, imag_sign) - 0.5) * 2; 
            dist = the_dist_sign .* abs(Rotated_fiber_CoM_xy_1st - PAs_xy);
            % Find positive elements
            positive_indices = find(dist > 0);
            % Find the minimum positive value
            min_positive_value = min(dist(positive_indices));
            % Find the position of the minimum positive value
            min_positive_position = find(dist == min_positive_value);
            ref_PAs_xy = PAs_xy(min_positive_position);

            Rotated_fiber_CoM_xy_normalized = [Rotated_fiber_CoM_xy(1,:) - real(ref_PAs_xy); ...
                Rotated_fiber_CoM_xy(2,:) - imag(ref_PAs_xy)];

%             plot(Rotated_fiber_CoM_xy_normalized(1,:), Rotated_fiber_CoM_xy_normalized(2,:));; hold on

            save([save_pathname, filesep, 'PlusInfo_', filename(14: end-4), '.mat'],'centers', ...
                'circleMask', 'CoM_x', 'Energy', 'Good_case_frm', 'L_ee_norm', ...
                'L_ee_norm_belowOne', 'lzero','metric','radii','xy', 'Chi', 'aniso', ...
                'Rotated_fiber_CoM_xy_normalized', 'RotMatrix_correct', 'Ctr2Ctr_x', 'Ctr2Ctr_y');

            clearvars CoM_x Energy xy centers radii L_ee_norm_belowOne L_ee_norm Chi aniso Rotated_fiber_CoM_xy_normalized
        end
    end
end
