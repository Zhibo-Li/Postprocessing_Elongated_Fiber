%%%%%%%%%%%%%%%%%% calculate contact forces %%%%%%%%%%%%%%%%%%%%%%%%%

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
for sub1Path_i = 7
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

            ster_forces_y = zeros(length(fileinfo), 1);
            for ii = 1:length(fileinfo)

                currentVTK_name = fileinfo(ii).name;
                steric_fib_obs = strcat('fiber_beads_f_steric_fib_obs_', currentVTK_name(end-9: end-4), '.dat');

                if str2double(currentVTK_name(end-9: end-4)) ~= 0

                    ster_forces = readmatrix(fullfile(fileinfo(ii).folder, steric_fib_obs));
%                     abs_ster_forces = sqrt(sum(ster_forces.^2, 2));
                    ster_forces_y(ii) = -sum(ster_forces(:, 1));
                end
            end

            if sum(ster_forces_y) ~= 0
                ster_forces_y(ster_forces_y == 0) = [];
                impulse_time_integral = mean(ster_forces_y)*numel(ster_forces_y);
                mean_force = mean(ster_forces_y);
            else
                impulse_time_integral = 0;
                mean_force = 0;
            end

            tmp_index =  thetheta0==deg_num & theL==L_num & ismembertol(they0, y0_num,0.0125);
            Ind = find(tmp_index==1);

            Loc = ['L', num2str(Ind+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
            writematrix(impulse_time_integral,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value into the excel.

%             Loc = ['K', num2str(Ind+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
%             writematrix(mean_force,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value into the excel.

            clearvars impulse_time_integral  

        end
    end
end


%% plot
clear; close all; clc;

xlsfile = readcell(['D:\Dropbox\Collaboration - LadHyX\' ...
    'Give_to_Zhibo_nonShared\Data_Give_to_Zhibo_20230223\' ...
    'results_2023_02_23_test.xlsx'],'Sheet','Sheet1','NumHeaderLines',1); 

thedyn = cell2mat(xlsfile(:, 5));

theL = cell2mat(xlsfile(:, 2)); %theL(thedyn==3) = [];
thedelta = cell2mat(xlsfile(:, 4)); thedelta(thedyn==3) = [];
theimpulse = cell2mat(xlsfile(:, 9)); %theimpulse(thedyn==3) = [];
theforce = cell2mat(xlsfile(:, 11)); %theforce(thedyn==3) = [];
theforce_times_t = cell2mat(xlsfile(:, 12)); theforce_times_t(thedyn==3) = [];


% figure; plot(theforce, 'ro')
% figure; plot(thedelta, 'b*')
figure; plot(theforce_times_t, thedelta, 'm^')
% xlim([-0.5e-4 0.5e-4])

% xlabel('$Steric\ force\ (f_y)$','FontSize', 22,'Interpreter', 'latex'); 
xlabel('$\int f_y t$','FontSize', 22,'Interpreter', 'latex'); 
% xlabel('$Interaction\ index$','FontSize', 22,'Interpreter', 'latex'); 
ylabel('$Deviation\ (\delta/h_{\rm obs})$','FontSize', 22,'Interpreter', 'latex');
set_plot(gcf, gca)
