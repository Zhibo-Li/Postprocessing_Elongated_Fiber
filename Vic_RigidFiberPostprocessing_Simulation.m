clear; close all; clc;

[excelname, excelpathname] = uigetfile(['D:\Dropbox\Collaboration - LadHyX\' ...
    'Give_to_Zhibo_nonShared\Data_Give_to_Zhibo_20230223\.xlsx'], ...
    'Please choose the excel accordingly and carefully!!');
xlsfile = readcell([excelpathname, excelname],'Sheet','Sheet1','NumHeaderLines',1); 
thedeg = cell2mat(xlsfile(:, 1)); 
theL = round(cell2mat(xlsfile(:, 2)), 1); 
they0 = round(cell2mat(xlsfile(:, 3)), 3); 

% obs_beads_size = 0.6e-6; % notice here (radius)
fiber_beads_size = 2e-6; % notice here (radius)

load(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\input_data\obs_2d_20230223.mat']); 
% get this file based on the following code in 'plot contact infromation vs
% initial condition' part.
load(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\input_data\pert_cont (50%).mat'])
pert_C1_50 = pert_C1; pert_C2_50 = pert_C2;
load(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\input_data\pert_cont (30%).mat'])
pert_C1_30 = pert_C1; pert_C2_30 = pert_C2;
l_obstacle = 86.6; % um
time_step = 0.01; % s

parent_path = 'D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\Data_Give_to_Zhibo_20230223\simulations';
sub1_path = dir(parent_path);
for sub1Path_i = 3:length(sub1_path)
    current_deg = sub1_path(sub1Path_i).name; 
    newStr = strrep(current_deg,'o','.'); newStr = strrep(newStr,'m','-'); 
    deg_num = str2double(newStr(8:end));
    sub2_path = dir(fullfile(parent_path, current_deg));
    for sub2Path_i = 3:length(sub2_path)
        current_L = sub2_path(sub2Path_i).name;
        newStr = strrep(current_L,'o','.');
        L_num = str2double(newStr(3:end));
        sub3_path = dir(fullfile(parent_path, current_deg, current_L));
        for sub3Path_i = 3:length(sub3_path)
            current_y0 = sub3_path(sub3Path_i).name;
            newStr = strrep(current_y0,'o','.');
            y0_num = str2double(newStr(4:end));
            if isnan(y0_num)
                y0_num = str2double(newStr(4:end-10));
            end
            fileinfo = dir(fullfile(parent_path, current_deg, current_L, current_y0, 'output_data\*.vtk'));

            perturb_contact_50 = 0; perturb_contact_30 = 0; direct_contact = 0;
            for ii = 1:length(fileinfo)
                snapshot = readVTK(fullfile(fileinfo(ii).folder, fileinfo(ii).name));
                XY = snapshot.points(:, 1:2);

                if_direct_contact = min(pdist2(XY,obs_2d,'euclidean','Smallest',1)) < fiber_beads_size * 1.1; % if it's direct contact.
                if if_direct_contact
                    direct_contact = direct_contact + 1;
                end

                %%%%%% 50%
                in_triangle_50 = inpolygon(XY(:,1),XY(:,2),obs_2d(:,1),obs_2d(:,2)); % the points on the fiber which are inside the triangular obstacle
                in_perturb_C1_50 = inpolygon(XY(:,1),XY(:,2),pert_C1_50(:,1),pert_C1_50(:,2)); % the points on the fiber which are inside the perturbed line 1
                if_inbetween_twolines_50 = and(~logical(in_triangle_50), logical(in_perturb_C1_50)); % the fiber is in between the two contours
                if_inbetween_twolines_50 = logical(sum(if_inbetween_twolines_50));

                if min(pdist2(XY,pert_C1_50,'euclidean','Smallest',1)) < fiber_beads_size * 1.1 || ...
                    if_direct_contact || if_inbetween_twolines_50
                    perturb_contact_50 = perturb_contact_50 + 1;
                end

                in_perturb_C2_50 = inpolygon(XY(:,1),XY(:,2),pert_C2_50(:,1),pert_C2_50(:,2)); % the points on the fiber which are inside the perturbed line 2
                if_in_C2_50 = logical(sum(logical(in_perturb_C2_50)));

                if (min(pdist2(XY,pert_C2_50,'euclidean','Smallest',1)) < fiber_beads_size * 1.1 || if_in_C2_50) && ~if_direct_contact % counting number of perturbed contact
                    perturb_contact_50 = perturb_contact_50 + 1;
                end

                %%%%%% 30%
                in_triangle_30 = inpolygon(XY(:,1),XY(:,2),obs_2d(:,1),obs_2d(:,2)); % the points on the fiber which are inside the triangular obstacle
                in_perturb_C1_30 = inpolygon(XY(:,1),XY(:,2),pert_C1_30(:,1),pert_C1_30(:,2)); % the points on the fiber which are inside the perturbed line 1
                if_inbetween_twolines_30 = and(~logical(in_triangle_30), logical(in_perturb_C1_30)); % the fiber is in between the two contours
                if_inbetween_twolines_30 = logical(sum(if_inbetween_twolines_30));

                if min(pdist2(XY,pert_C1_30,'euclidean','Smallest',1)) < fiber_beads_size * 1.1 || ...
                    if_direct_contact || if_inbetween_twolines_30
                    perturb_contact_30 = perturb_contact_30 + 1;
                end

                in_perturb_C2_30 = inpolygon(XY(:,1),XY(:,2),pert_C2_30(:,1),pert_C2_30(:,2)); % the points on the fiber which are inside the perturbed line 2
                if_in_C2_30 = logical(sum(logical(in_perturb_C2_30)));

                if (min(pdist2(XY,pert_C2_30,'euclidean','Smallest',1)) < fiber_beads_size * 1.1 || if_in_C2_30) && ~if_direct_contact % counting number of perturbed contact
                    perturb_contact_30 = perturb_contact_30 + 1;
                end

            end

            interaction1 = direct_contact * time_step / (l_obstacle / U_max_phy);
            interaction2 = perturb_contact_50 * time_step / (l_obstacle / U_max_phy);
            interaction3 = perturb_contact_30 * time_step / (l_obstacle / U_max_phy);

            tmp_index =  thedeg==deg_num & theL==L_num & ismembertol(they0, y0_num,0.0125);
            Ind = find(tmp_index==1);
            Loc = ['M', num2str(Ind+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
            writematrix(interaction1,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value inti the excel.
            Loc = ['N', num2str(Ind+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
            writematrix(interaction2,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value inti the excel.
            Loc = ['O', num2str(Ind+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
            writematrix(interaction3,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value inti the excel.
            
        end
    end
end





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% calculate contact infromation %%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

[excelname, excelpathname] = uigetfile(['D:\Dropbox\Collaboration - LadHyX\' ...
    'Give_to_Zhibo_nonShared\Data_Give_to_Zhibo_20230223\.xlsx'], ...
    'Please choose the excel accordingly and carefully!!');
xlsfile = readcell([excelpathname, excelname],'Sheet','Sheet1','NumHeaderLines',1); 
thetheta0 = cell2mat(xlsfile(:, 1)); 
theL = round(cell2mat(xlsfile(:, 2)), 1); 
they0 = round(cell2mat(xlsfile(:, 3)), 8); 

obs_beads_size = 0.6e-6; % notice here (radius)
fiber_beads_size = 2e-6; % notice here (radius)

% %%% obstacle position and save it as *.mat file
obs = readVTK('D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\Data_Give_to_Zhibo_20230223\input_data\obstacle_beads.vtk');
% obs_beads_size = 0.6e-6; % notice here (radius)
obs_2d = obs.points(1:208, 1:2); obs_2d_center = mean(obs_2d, 1);  % 208 because there are multiple layers.
% sort the coordinates clockwise
[theta, ~] = cart2pol(obs_2d(:,1)-obs_2d_center(1), obs_2d(:,2)-obs_2d_center(2));
obs_2d = sortrows([obs_2d, theta], 3); obs_2d = obs_2d(:, 1:2);

ver_ind = find(obs_2d(:,2) == max(obs_2d(:,2))); % index of the vertex
obs_2d = [obs_2d(ver_ind:end, :); obs_2d(1:ver_ind-1, :)]; % re-order the coordinates (highest apex first).
windward_edge_x = obs_2d(1:round(size(obs_2d, 1)/3)+3,1); % 1/3 of all the points; +3 is to extend the edge a bit more.
windward_edge_y = obs_2d(1:round(size(obs_2d, 1)/3)+3,2);
% the highest and lowest of the obstacle for y_c calculation.
obs_top = max(obs_2d(:,2));
obs_bottom = min(obs_2d(:,2));

parent_path = 'D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\Data_Give_to_Zhibo_20230223\simulations';
sub1_path = dir(parent_path);
for sub1Path_i = 3:length(sub1_path)
    current_deg = sub1_path(sub1Path_i).name; 
    newStr = strrep(current_deg,'o','.'); newStr = strrep(newStr,'m','-'); 
    deg_num = str2double(newStr(8:end));
    sub2_path = dir(fullfile(parent_path, current_deg));
    for sub2Path_i = 3:length(sub2_path)
        current_L = sub2_path(sub2Path_i).name;
        newStr = strrep(current_L,'o','.');
        L_num = str2double(newStr(3:end));
        sub3_path = dir(fullfile(parent_path, current_deg, current_L));
        for sub3Path_i = 3:length(sub3_path)
            current_y0 = sub3_path(sub3Path_i).name;
            newStr = strrep(current_y0,'o','.');
            y0_num = str2double(newStr(4:end));
            if isnan(y0_num)
                y0_num = str2double(newStr(4:end-10));
            end
            fileinfo = dir(fullfile(parent_path, current_deg, current_L, current_y0, 'output_data\*.vtk'));

            for ii = 1:length(fileinfo)
                snapshot = readVTK(fullfile(fileinfo(ii).folder, fileinfo(ii).name));
                XY = snapshot.points(:, 1:2);

                the_dist = pdist2(XY,[windward_edge_x, windward_edge_y],'euclidean');

                if min(min(the_dist)) < 0.593 * (obs_beads_size+fiber_beads_size)  % if it's direct contact.

                    com_XY = mean(XY, 1); % CoM of current fiber
                    Gyr = 1/size(XY,1) * [sum((XY(:, 1)-com_XY(1)).^2),  sum((XY(:, 1)-com_XY(1)) .* (XY(:, 2)-com_XY(2)));
                        sum((XY(:, 2)-com_XY(2)) .* (XY(:, 1)-com_XY(1))), sum((XY(:, 2)-com_XY(2)).^2)];
                    [eigenV,eigenD] = eig(Gyr);
                    [d,ind] = sort(diag(eigenD));
                    Ds = eigenD(ind,ind);
                    Vs = eigenV(:,ind);
                    theta_c = atan(Vs(2,2)/Vs(1,2))/pi*180; % calculate theta_c

                    [minDist_ind_fiber, minDist_indy_obs] = find(the_dist == min(min(the_dist))); 
                    y_c = (windward_edge_y(minDist_indy_obs)-obs_bottom) / (obs_top-obs_bottom); % calculate y_c

                    if minDist_ind_fiber == 1 || minDist_ind_fiber == size(XY, 1)
                        if_fiberENDs_contact = 1; % check if the fiber ends contact the obstacle
                    else
                        if_fiberENDs_contact = 0;
                    end

                    ite_contact = ii - 1;  % the frame number of the contact case

                    break
                end
            end

            if exist('theta_c','var')
                tmp_index =  thetheta0==deg_num & theL==L_num & ismembertol(they0, y0_num,0.0125);
                Ind = find(tmp_index==1);

                Loc = ['I', num2str(Ind+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
                writematrix(y_c,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value into the excel.
                Loc = ['J', num2str(Ind+1)];  
                writematrix(theta_c,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  

                Loc = ['K', num2str(Ind+1)]; 
                writematrix(ite_contact,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  
                % load orientation information which is calculated from: Vic_RigidFiberPostprocessing_Simulation_orientations.m
                load(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared' ...
                    '\Data_Give_to_Zhibo_20230223_videos\Orientations\', current_deg, '_', current_L, ...
                    '_', current_y0, '_orientations.mat']);
                Loc = ['L', num2str(Ind+1)]; 
                writematrix(ori_ee(ite_contact),[excelpathname, excelname],'Sheet','Sheet1','Range', Loc); 

                Loc = ['P', num2str(Ind+1)]; 
                writematrix(if_fiberENDs_contact,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  
            end

            clearvars theta_c y_c ite_contact if_fiberENDs_contact

        end
    end
end


