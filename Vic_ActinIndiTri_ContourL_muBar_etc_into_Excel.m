%%%% Calculate contour length and mu_bar (actin vs. triangular pillar)
%
%
% data from Vic_ActinIndiTri_BendingE_Lee_Chi_Aniso.m
% data name format: PlusInfo_trajectory_..._.mat
% results save in the Excels named: xxxxxxxx-Actin-Individual_triangularPillar_uppoint.xlsx

%%%% The contour L is the average of 10% longest snapshots except for the
%%%% extreme values (out of 2-sigma among those averaged ones). 

clear; close all; clc;

mag = 0.1; % um/pixel

% load the data
data_path = uigetdir('F:\Processing & Results\FSI - Actin &  Individual Obstacle\', ...
    'Choose the data path which contains the *.mat results !');
Files = dir([data_path,'\*.mat']);

[excelname, excelpathname] = uigetfile(['F:\Processing & Results\FSI - Actin &  Individual Obstacle' ...
    '\Exp data triangular pillar lateral pointing\*.xlsx'], 'Choose right Excel !');
xlsfile_to_write = readcell([excelpathname, excelname],'Sheet','Sheet1','NumHeaderLines',1);

% load slope correction data to correct delta
correct_ind = find(cellfun(@(x) contains(x, 'PlusInfo_Slope_correction'), {Files.name}));
load(fullfile(Files(1).folder, Files(correct_ind).name));
% calculate the information of the obstacle
obs_2d(:, 2) = 2048 - obs_2d(:, 2); % !!!!!! Now, it's in matrix coordinates and the triangle pointing down !!!.
% !!!!!! and the filament data is also in matrix coordinates !
obs_center_xy = mean(obs_2d);
D = sort(pdist2(obs_center_xy, obs_2d));  % distances from center to edge (sorted);
obs_h = mean(D(1:round(0.01*length(D)))) + mean(D(end-round(0.01*length(D)):end));
% obstacle altitude is sum of 1% shortest average and 1% longest average values.
obs_base_y = obs_center_xy(2) + mean(D(1:round(0.01*length(D))));
% obstacle base_y is center_y plus 1% shortest average.

left_border = obs_center_xy(1) - 3.5*obs_h;
right_border = obs_center_xy(1) + 3.5*obs_h; % give range for y_0 and delta calculation
% correction value of delta
delta_correct = mean(CoM_y(CoM_x > right_border)) - mean(CoM_y(CoM_x < left_border));

all_names = xlsfile_to_write(:, 1);

for file_ind = 1:length(Files)

    filename = Files(file_ind).name;

    if contains(filename, 'PlusInfo_trajectory')

        load(fullfile(Files(1).folder, Files(file_ind).name));

%         figure; plot(CoM_x, CoM_y, 'ok'); hold on; plot(obs_2d(:,1), obs_2d(:,2), 'r'); 
%         axis equal

        excel_pos = find(cellfun(@(x) contains(filename(10:end-4), x), all_names));

        ContourL_all = xy.arclen_spl(Good_case_frm);
        ContourL = VicFc_Get_ContourLength(ContourL_all) * mag; % unit: um

        try
            y_0 = mean(CoM_y(CoM_x < left_border)); y_f = mean(CoM_y(CoM_x > right_border));
            theta_0 = Chi(1);
        catch
            y_0 = nan; y_f = nan; theta_0 = nan;
        end

        delta = -((y_f - y_0) - delta_correct); % delta is respect to the triangle pointing up.
        
        Loc = ['B', num2str(excel_pos+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
        writematrix(ContourL,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value inti the excel.
        Loc = ['C', num2str(excel_pos+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
        writematrix((obs_base_y - y_0)/obs_h,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value inti the excel.
        Loc = ['D', num2str(excel_pos+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
        writematrix(delta/obs_h,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value inti the excel.
        Loc = ['E', num2str(excel_pos+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
        writematrix(-atand(theta_0),[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value inti the excel.
        
    end
end
