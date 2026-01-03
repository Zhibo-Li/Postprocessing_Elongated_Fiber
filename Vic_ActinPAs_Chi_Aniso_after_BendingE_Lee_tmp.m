%%%% Calculate the orientation and sphericity of actin filaments and 
%%%% and add to the results after BendingE and Lee.
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

            Chi = zeros(size(Good_case_frm,2), 1);
            aniso =  zeros(size(Good_case_frm,2), 1);
            for frm_ind = 1:size(Good_case_frm,2)

                xy_ind = Good_case_frm(frm_ind);% index of the 'good' cases

                crd = xy.crd{1,xy_ind};
                CoM_xy = xy.centroid{1,xy_ind};

                Gyr = 1/size(crd,1) * [sum((crd(:, 1)-CoM_xy(1)).^2),  sum((crd(:, 1)-CoM_xy(1)) .* (crd(:, 2)-CoM_xy(2)));
                    sum((crd(:, 2)-CoM_xy(2)) .* (crd(:, 1)-CoM_xy(1))), sum((crd(:, 2)-CoM_xy(2)).^2)];

                [eigenV,eigenD] = eig(Gyr);
                [d,ind] = sort(diag(eigenD));
                Ds = eigenD(ind,ind);
                Vs = eigenV(:,ind);
                Chi(frm_ind) = atan(Vs(2,2)/Vs(1,2));

                Lambda1 = eigenD(2,2); Lambda2 =  eigenD(1,1);
                aniso(frm_ind) = 1 - 4*Lambda1*Lambda2/(Lambda1+Lambda2)^2;
              
            end

            save([save_pathname, filesep, 'PlusInfo_', filename(14: end-4), '.mat'],'centers', ...
                'circleMask', 'CoM_x', 'Energy', 'Good_case_frm', 'L_ee_norm', ...
                'L_ee_norm_belowOne', 'lzero','metric','radii','xy', 'Chi', 'aniso');

            clearvars CoM_x Energy xy centers radii L_ee_norm_belowOne L_ee_norm Chi aniso 
        end
    end
end
