%%%% Add information to the trajectory_M63_..._batch1.mat of filament.
% 1. Good_case_frm
% 2. obstacle locations
%%%% because the information is different from one case to another.

clear; close all; clc

%%%
scalefactor = 0.1;  % um/pix
delta_t = 0.02; % second

% load the obstacle position
[obs_file, obs_path] = uigetfile(['F:\Processing & Results\' ...
    'FSI - Actin &  Individual Obstacle\*.mat'], ...
    'Choose the corresponding obstacle !');

load(fullfile(obs_path, obs_file));
[row,col] = ind2sub(size(tri),find(tri == max(max(tri)))); 
obs_2d = [col, row]; obs_2d_center = mean(obs_2d, 1); 

% sort the coordinates clockwise
[theta, ~] = cart2pol(obs_2d(:,1)-obs_2d_center(1), obs_2d(:,2)-obs_2d_center(2));
obs_2d = sortrows([obs_2d, theta], 3); obs_2d = obs_2d(:, 1:2);

ver_ind = find(obs_2d(:,2) == max(obs_2d(:,2))); % index of the vertex
obs_2d = [obs_2d(ver_ind:end, :); obs_2d(1:ver_ind-1, :)]; % re-order the coordinates (highest apex first)
%%% should notice that the positions are in the 'image' coordinates.

% load the data
data_path = uigetdir('F:\Processing & Results\FSI - Actin &  Individual Obstacle\', ...
    'Choose the data path which contains the *.mat results !');
Files = dir([data_path,'\*.mat']);

for case_No = 1: length(Files)

    load(fullfile(Files(1).folder, Files(case_No).name))
    Good_case_frm = find(ismember(xy(1).frame, Good_case));

    save_filename = ['AddInfo_', Files(case_No).name];

    %%% save ...
    save([data_path, filesep , save_filename], 'framelist', ...
        'Good_case','InfoImage','prcs_img','prmt','ROI','xy', ...
        'Good_case_frm','obs_2d')

end



