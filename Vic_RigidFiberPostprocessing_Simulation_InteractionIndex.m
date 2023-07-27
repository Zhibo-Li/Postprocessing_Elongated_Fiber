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

            Loc = ['K', num2str(Ind+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
            writematrix(mean_force,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value into the excel.

            % interaction index 3 (defines as contact probability -- varying layer thickness)
            circle_intersec = 0; line_intersec = 0;
            for jj = 1:length(fileinfo)

                snapshot = readVTK(fullfile(fileinfo(jj).folder, fileinfo(jj).name));
                XY = snapshot.points(:, 1:2);
                Half_L = sqrt((XY(1,1)-XY(end,1))^2 + (XY(1,2)-XY(end,2))^2) / 2 + 0.593 * (obs_beads_size+fiber_beads_size);

                theta = linspace(0,2*pi,300);
                circle_x = Half_L*cos(theta); circle_x = circle_x + mean(XY(:, 1));
                circle_y = Half_L*sin(theta); circle_y = circle_y + mean(XY(:, 2));
                [intersec_cir_x, intersec_cir_y] = polyxpoly(obs_2d(:,1), obs_2d(:,2), circle_x, circle_y);
                if ~isempty(intersec_cir_y)
                    circle_intersec = circle_intersec + 1;
                    if min(pdist2(XY,obs_2d,'euclidean','Smallest',1)) < 0.593 * (obs_beads_size+fiber_beads_size)
                        line_intersec = line_intersec + 1;
                    end
                end
            end
            if circle_intersec == 0
                interaction3 = 0;
            else
                interaction3 = line_intersec / circle_intersec;
            end

            tmp_index =  thetheta0==deg_num & theL==L_num & ismembertol(they0, y0_num,0.0125);
            Ind = find(tmp_index==1);

            Loc = ['M', num2str(Ind+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
            writematrix(interaction3 ,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value into the excel

            clearvars interaction3

        end
    end
end


%% plot steric force vs delta
clear; close all; clc;

xlsfile = readcell(['D:\Dropbox\Collaboration - LadHyX\' ...
    'Give_to_Zhibo_nonShared\Data_Give_to_Zhibo_20230223\' ...
    'results_2023_02_23_test.xlsx'],'Sheet','Sheet1','NumHeaderLines',1); 

thedyn = cell2mat(xlsfile(:, 5)); % the dynamics: 3 is trapping

thedelta = cell2mat(xlsfile(:, 4)); % thedelta(thedyn==3) = [];
theforce = cell2mat(xlsfile(:, 11)); % theforce(thedyn==3) = [];

figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
plot(theforce / max(abs(theforce)), thedelta, '.k', 'MarkerSize', 15)

xlim([-0.1 1.1]); ylim([-0.3 0.7]);
xlabel('$f_y/|f_y|_{\rm max}$','FontSize', 24,'Interpreter', 'latex'); 
ylabel('$\delta$','FontSize', 24,'Interpreter', 'latex');
text(-0.05, 0.6, 'Simulation','FontSize', 24, 'Interpreter', 'latex','BackgroundColor',[.7 .7 .7])
set(gca,'Box', 'On','XGrid', 'On','YGrid', 'On','FontSize', 24,'TickLabelInterpreter','latex')

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['F:\Processing & Results\FSI - Rigid Fiber &  Individual ' ...
    'Obstacle\Figures\about interaction index\F_y_delta_simu.pdf']);

%% plot interaction index vs delta
clear; close all; clc;

xlsfile = readcell(['D:\Dropbox\Collaboration - LadHyX\' ...
    'Give_to_Zhibo_nonShared\Data_Give_to_Zhibo_20230223\' ...
    'results_2023_02_23_test.xlsx'],'Sheet','Sheet1','NumHeaderLines',1); 

thedyn = cell2mat(xlsfile(:, 5)); % the dynamics: 3 is trapping

thedelta = cell2mat(xlsfile(:, 4)); % thedelta(thedyn==3) = [];
contact_prob = cell2mat(xlsfile(:, 13)); % theforce(thedyn==3) = [];

figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
plot(contact_prob, thedelta, '.k', 'MarkerSize', 15)

xlim([-0.1 1.1]); ylim([-0.3 0.7]);
xlabel('Contact probability','FontSize', 24,'FontName', 'Times New Roman'); 
ylabel('$\delta$','FontSize', 24,'Interpreter', 'latex');
text(-0.05, 0.6, 'Simulation','FontSize', 24, 'Interpreter', 'latex','BackgroundColor',[.7 .7 .7])
set(gca,'Box', 'On','XGrid', 'On','YGrid', 'On','FontSize', 24,'TickLabelInterpreter','latex')

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['F:\Processing & Results\FSI - Rigid Fiber &  Individual ' ...
    'Obstacle\Figures\about interaction index\Contact_probability_delta_simu.pdf']);

