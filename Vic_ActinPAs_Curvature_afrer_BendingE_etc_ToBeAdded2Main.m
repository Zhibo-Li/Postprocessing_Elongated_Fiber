%%%% Calculate the curvature and instantaneous speed of actin filaments 
%%%% then add to the results after BendingE, Lee, orientation and sphericity.
%
% data from PlusInfo_trajectory_..._batch1.mat
% saving name format: PlusInfo_trajectory_..._batch1.mat
%

clear; close all; clc;

mag = 0.1; % um/pixel

xlsfile = readcell('ForActinPostprocessing.xlsx','Sheet','Sheet1','NumHeaderLines',1);
% This is the file that contains all the information about the later processing (in sheet 1).

NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.
Init_U = xlsfile(:, 8);  % initial velocity, unit: m/s.

for no_Group = [7 8 13:28]

    the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
    thefiles = dir(fullfile(storePath{no_Group},'*.mat'));

    for file_ind = 1:length(thefiles)

        filename = thefiles(file_ind).name;

        if contains(filename, 'PAsInfoAdded_')

            pathname = thefiles(1).folder;
            save_pathname = strrep(pathname,'results','results_plus');
            filename = thefiles(file_ind).name
            load([save_pathname, filesep, 'PlusInfo_', filename(14: end-4), '.mat'])

            for frm_ind = 1:size(Good_case_frm,2)

                xy_ind = Good_case_frm(frm_ind);% index of the 'good' cases

                spl = xy.spl{1,xy_ind};
                [L2,R2,K2] = curvature(spl);
                R2(isnan(R2)) = inf;
                RR2 = movmean(R2,ceil(size(spl,1)/100));
                Curvature = 1./RR2;
%                 figure; scatter(spl(:,1), spl(:, 2),5, 1./RR2,"filled"); axis equal
%                 close
                xy.curvature{1,xy_ind} = Curvature; % unit: 1/pixel

            end

            save([save_pathname, filesep, 'PlusInfo_', filename(14: end-4), '.mat'], 'centers', ...
                'circleMask', 'CoM_x', 'CoM_y','Energy', 'Good_case_frm', 'L_ee_norm', ...
                'L_ee_norm_belowOne', 'lzero','metric','radii','xy', 'Chi', 'aniso', ...
                'Ctr2Ctr_x', 'Ctr2Ctr_y', 'Rotated_fiber_CoM_xy_normalized', ...
                'RotMatrix_correct', 'fiber_CoM_xy_in_lattice', 'Instan_V');

            clearvars xy CoM_x CoM_y Good_case_frm Instan_V Curvature
        end
    end
end
