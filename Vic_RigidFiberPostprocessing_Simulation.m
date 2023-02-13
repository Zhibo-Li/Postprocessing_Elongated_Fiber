clear; close all; clc;

[excelname, excelpathname] = uigetfile(['D:\Dropbox\Collaboration - LadHyX\' ...
    'Give_to_Zhibo_nonShared\.xlsx'], 'Please choose the excel accordingly and carefully!!');
xlsfile = readcell([excelpathname, excelname],'Sheet','Sheet1','NumHeaderLines',1); 
thedeg = cell2mat(xlsfile(:, 1)); 
theL = round(cell2mat(xlsfile(:, 2)), 1); 
they0 = round(cell2mat(xlsfile(:, 3)), 3); 

% obs = readVTK('D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\obstacle_beads.vtk');
% obs_beads_size = 1.5e-6; % notice here (radius)
fiber_beads_size = 0.8e-6; % notice here (radius)
% obs_2d = obs.points(1:267, 1:2); obs_2d_center = mean(obs_2d, 1);
% obs_2d = obs_2d * 8e-6;
% obs_2d_center = obs_2d_center * 8e-6;
% obs_2d(:, 1) = (obs_2d(:, 1) - obs_2d_center(1)) * ...
%     (2.5e-5 + obs_beads_size) / 2.5e-5 + obs_2d_center(1);
% obs_2d(:, 2) = (obs_2d(:, 2) - obs_2d_center(2)) * ...
%     (2.5e-5 + obs_beads_size) / 2.5e-5 + obs_2d_center(2); 
% %%% 2.5e-5 is the 1/3 the height of the equilateral triangle
% obs_2d = [obs_2d; obs_2d(1, :)];
load(['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\' ...
    'Figures\about interaction index\obs_2d.mat']);
load(['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\' ...
    'Figures\about interaction index\pert_cont (50%).mat'])
l_obstacle = 86.6; % um
time_step = 0.01; % s

parent_path = 'D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\data';
sub1_path = dir(parent_path);
for sub1Path_i = 3:length(sub1_path)
    current_deg = sub1_path(sub1Path_i).name; 
    newStr = strrep(current_deg,'o','.'); newStr = strrep(newStr,'m','-'); 
    deg_num = str2double(newStr(6:end));
    sub2_path = dir(fullfile(parent_path, current_deg));
    for sub2Path_i = 3:length(sub2_path)
        current_L = sub2_path(sub2Path_i).name;
        newStr = strrep(current_L,'o','.');
        L_num = str2double(newStr(13:end));
        sub3_path = dir(fullfile(parent_path, current_deg, current_L));
        for sub3Path_i = 3:length(sub3_path)
            current_y0 = sub3_path(sub3Path_i).name;
            newStr = strrep(current_y0,'o','.');
            y0_num = str2double(newStr(4:end));
            fileinfo = dir(fullfile(parent_path, current_deg, current_L, current_y0, 'output_data\*.vtk'));

            perturb_contact = 0; direct_contact = 0;
            for ii = 1:length(fileinfo)
                snapshot = readVTK(fullfile(fileinfo(ii).folder, fileinfo(ii).name));
                XY = snapshot.points(:, 1:2);

                if_direct_contect = min(pdist2(XY,obs_2d,'euclidean','Smallest',1)) < fiber_beads_size * 1.25; % if it's direct contact.
                if if_direct_contect
                    direct_contact = direct_contact + 1;
                end

%                 in_triangle = inpolygon(XY(:,1),XY(:,2),obs_2d(:,1),obs_2d(:,2)); % the points on the fiber which are inside the triangular obstacle
%                 in_perturb_C1 = inpolygon(XY(:,1),XY(:,2),pert_C1(:,1),pert_C1(:,2)); % the points on the fiber which are inside the perturbed line 1
%                 if_inbetween_twolines = and(~logical(in_triangle), logical(in_perturb_C1)); % the fiber is in between the two contours
%                 if_inbetween_twolines = logical(sum(if_inbetween_twolines));
% 
%                 if min(pdist2(XY,pert_C1,'euclidean','Smallest',1)) < fiber_beads_size * 1.25 || ...
%                     if_direct_contect || if_inbetween_twolines
%                     perturb_contact = perturb_contact + 1;
%                 end

                in_perturb_C2 = inpolygon(XY(:,1),XY(:,2),pert_C2(:,1),pert_C2(:,2)); % the points on the fiber which are inside the perturbed line 2
                if_in_C2 = logical(sum(logical(in_perturb_C2)));

                if (min(pdist2(XY,pert_C2,'euclidean','Smallest',1)) < fiber_beads_size * 1.25 || if_in_C2) && ~if_direct_contect % counting number of perturbed contact
                    perturb_contact = perturb_contact + 1;
                end

            end

            interaction1 = perturb_contact * time_step / (l_obstacle / U_max_phy);
            interaction2 = direct_contact * time_step / (l_obstacle / U_max_phy);
%             interaction3 = interaction1 + interaction2;

            tmp_index =  thedeg==deg_num & theL==L_num & they0==y0_num;
            Ind = find(tmp_index==1);
            Loc = ['AD', num2str(Ind+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
            writematrix(interaction1,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value inti the excel.
            Loc = ['AE', num2str(Ind+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
            writematrix(interaction2,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value inti the excel.
%             Loc = ['X', num2str(Ind+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
%             writematrix(interaction3,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value inti the excel.
            
        end
    end
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
xlsfile = readcell(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'results_more-data.xlsx'],'Sheet','Sheet1','NumHeaderLines',1); 
together_plot = [cell2mat(xlsfile(:, 1:4)), cell2mat(xlsfile(:, 13)), cell2mat(xlsfile(:, 27))]; 
together_plot = together_plot';
% thedeg = cell2mat(xlsfile(:, 1)); 
% theL = round(cell2mat(xlsfile(:, 2)), 1); 
% they0 = round(cell2mat(xlsfile(:, 3)), 3); 
% delta = cell2mat(xlsfile(:, 4));
% interaction2 = cell2mat(xlsfile(:, 11)); 
prompt = {'The lower bound of the initial angle:', 'The upper bound of the initial angle:', ...
    'The lower bound of the contour length:','The upper bound of the contour length:'...
    'The lower bound of the initial position:','The upper bound of the initial position:'};
definput = {'-10', '10', '0.5', '1.5', '0', '1'};
answer = inputdlg(prompt, 'Input (please input NaN if there is no bound)', [1 35] , definput);

%%%%%%%%%%%%%%%%%%%%%%%%% assign the values %%%%%%%%%%%%%%%%%%%%%%%%%%
range_chi0_low = str2double(answer{1,1}); if isnan(range_chi0_low); range_chi0_low = -91; end
range_chi0_up = str2double(answer{2,1}); if isnan(range_chi0_up); range_chi0_up = 91; end
range_L_low = str2double(answer{3,1}); if isnan(range_L_low); range_L_low = 0; end
range_L_up = str2double(answer{4,1}); if isnan(range_L_up); range_L_up = 10; end
range_y0_low = str2double(answer{5,1}); if isnan(range_y0_low); range_y0_low = -10; end
range_y0_up = str2double(answer{6,1}); if isnan(range_y0_up); range_y0_up = 10; end
% chi_0
together_plot(:, together_plot(1, :) < range_chi0_low) = []; 
together_plot(:, together_plot(1, :) > range_chi0_up) = [];
% L
together_plot(:, together_plot(2, :) < range_L_low) = [];
together_plot(:, together_plot(2, :) > range_L_up) = [];
% y_0
together_plot(:, together_plot(3, :) < range_y0_low) = [];  
together_plot(:, together_plot(3, :) > range_y0_up) = [];


figure('color', 'w'); set(gcf, 'Position', [100 100 1000 600]);
% cmap = cmocean('thermal');
cmap = colormap("jet");

bypass_edge_together = together_plot(:, together_plot(5, :)==0); 
bypass_tip_together = together_plot(:, together_plot(5, :)==1);
pole_vaulting_together = together_plot(:, together_plot(5, :)==2); 
trapped_together = together_plot(:, together_plot(5, :)==3);

scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond'); hold on  % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'square'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', '^'); hold on % for legend only

scatter(trapped_together(6, :), trapped_together(4, :), 220, trapped_together(3, :), 'Filled', 'diamond','MarkerEdgeColor','k'); hold on 
scatter(bypass_edge_together(6, :), bypass_edge_together(4, :), 200, bypass_edge_together(3, :), 'Filled','MarkerEdgeColor','k'); hold on
scatter(bypass_tip_together(6, :), bypass_tip_together(4, :), 200, bypass_tip_together(3, :), 'Filled', 'square','MarkerEdgeColor','k'); hold on
scatter(pole_vaulting_together(6, :), pole_vaulting_together(4, :), 220, pole_vaulting_together(3, :), 'Filled', '^','MarkerEdgeColor','k'); hold on

% cmap(size(together_plot,2)); 
hcb=colorbar; caxis([0 1])
title(hcb,'$Initial\ position\ (y_0/h_{obs})$','FontSize', 16,'Interpreter', 'latex'); grid on
set(gca,'FontSize',16);
% xlabel('$max\left|U_0-U(t)\right|/U_0$','FontSize', 22,'Interpreter', 'latex');
% xlabel('$Contact\ probability\ (disturbed\ layer = 40\mu{m}) $','FontSize', 22,'Interpreter', 'latex');  % for interaction3
% xlabel('$Normalized\ direct\ contact\ duration$','FontSize', 22,'Interpreter', 'latex');  
xlabel('$Normalized\ perturbed\ duration$','FontSize', 22,'Interpreter', 'latex'); 
% xlabel('$Interaction\ index$','FontSize', 22,'Interpreter', 'latex'); 
ylabel('$Deviation\ (\delta/h_{obs})$','FontSize', 22,'Interpreter', 'latex');
legend({'Trapping','Below','Above','Pole-vaulting'}, 'Location', 'northeast','FontSize', 14,'Interpreter', 'latex')
% title('$0.5<L<1$','FontSize', 22,'Interpreter', 'latex')
% xlim([0 12]); %ylim([-0.4 0.8])
% f=gcf;
% exportgraphics(f,'F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\Figures\about interaction index\perturbed duration (0o30).png','Resolution',100)




%% ???
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% plot contact infromation vs initial condition %%%%%%%%%%%%%%%
%%% based on the information in excel from Clement (need to be checked)

clear; close all; clc;
xlsfile = readcell(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'results_2023_01_24.xlsx'],'Sheet','Sheet1','NumHeaderLines',1); 

%%%%%%%%%%%%%%%%%%%%%%%%% calculate the y_c %%%%%%%%%%%%%%%%%%%%%%%%%%
load(['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\' ...
    'Figures\about interaction index\obs_2d.mat']);

ver_ind = find(obs_2d(:,2) == max(obs_2d(:,2))); % index of the vertex
windward_edge_x = obs_2d(ver_ind:end,1); windward_edge_y = obs_2d(ver_ind:end,2);
% the edge of the obstacle to calculate y_c.

initial_info = [cell2mat(xlsfile(:, 1:3)), cell2mat(xlsfile(:, 5:6))]; 

for case_no = 1:size(initial_info, 1)

    theta0 = initial_info(case_no, 1);
    theta0_str = strrep(num2str(theta0),'.','o'); theta0_str = strrep(num2str(theta0_str),'-','m'); 
    L0 = initial_info(case_no, 2);
    L0_str = strrep(num2str(L0),'.','o');
    y0 = initial_info(case_no, 3);
    y0 = num2str(y0,'%5.3f');
    y0_str = strrep(num2str(y0),'.','o');
    theta_c = initial_info(case_no, 4);
    frame_no = initial_info(case_no, 5);
 
    filename =['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
        'simulations_latest\theta',theta0_str, '\L_over_l_obs', L0_str, ...
        '\y0_', y0_str, '\output_data\fibers_000', num2str(frame_no,'%03.f'), '.vtk'];

    snapshot = readVTK(filename);
    XY = snapshot.points(:, 1:2);

    L_current = sqrt((XY(1,1)-XY(end,1))^2 + (XY(1,2)-XY(end,2))^2);

    com_XY = mean(XY, 1); 
    L_start = [com_XY(1) - cosd(theta_c) * L_current, com_XY(2) - sind(theta_c)  * L_current];
    L_end = [com_XY(1) + cosd(theta_c) * L_current, com_XY(2) + sind(theta_c)  * L_current];
            
%     plot(windward_edge_x, windward_edge_y); hold on;
%     plot([L_start(1) L_end(1)], [L_start(2) L_end(2)]); hold on
%     plot(XY(:,1), XY(:,2))
%     close;

end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% plot contact infromation vs initial condition %%%%%%%%%%%%%%%
%%% based on simulation data (in 'simulations_latest' folder)
clear; close all; clc;

[excelname, excelpathname] = uigetfile(['D:\Dropbox\Collaboration - LadHyX\' ...
    'Give_to_Zhibo_nonShared\.xlsx'], 'Please choose the excel accordingly and carefully!!');
xlsfile = readcell([excelpathname, excelname],'Sheet','Sheet1','NumHeaderLines',1); 
thetheta0 = cell2mat(xlsfile(:, 1)); 
theL = round(cell2mat(xlsfile(:, 2)), 1); 
they0 = round(cell2mat(xlsfile(:, 3)), 8); 

% obs_beads_size = 1.5e-6; % notice here (radius)
fiber_beads_size = 0.8e-6; % notice here (radius)

load(['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\' ...
    'Figures\about interaction index\obs_2d.mat']);
ver_ind = find(obs_2d(:,2) == max(obs_2d(:,2))); % index of the vertex
windward_edge_x = obs_2d(ver_ind:end,1); windward_edge_y = obs_2d(ver_ind:end,2);
% the edge of the obstacle to calculate y_c.
L_edge = sqrt((windward_edge_x(1)-windward_edge_x(end))^2 + (windward_edge_y(1)-windward_edge_y(end))^2);
L_edge_X = mean(windward_edge_x);  L_edge_Y = mean(windward_edge_y);
L_start_obsedge = [L_edge_X - cosd(60) * L_edge, L_edge_Y - sind(60)  * L_edge];
L_end_obsedge = [L_edge_X + cosd(60) * L_edge, L_edge_Y + sind(60)  * L_edge];

parent_path = 'D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\simulations_latest';
sub1_path = dir(parent_path);
for sub1Path_i = 3:length(sub1_path)
    current_deg = sub1_path(sub1Path_i).name; 
    newStr = strrep(current_deg,'o','.'); newStr = strrep(newStr,'m','-'); 
    deg_num = str2double(newStr(6:end));
    sub2_path = dir(fullfile(parent_path, current_deg));
    for sub2Path_i = 3:length(sub2_path)
        current_L = sub2_path(sub2Path_i).name;
        newStr = strrep(current_L,'o','.');
        L_num = str2double(newStr(13:end));
        sub3_path = dir(fullfile(parent_path, current_deg, current_L));
        for sub3Path_i = 3:length(sub3_path)
            current_y0 = sub3_path(sub3Path_i).name;
            newStr = strrep(current_y0,'o','.');
            y0_num = str2double(newStr(4:end));
            fileinfo = dir(fullfile(parent_path, current_deg, current_L, current_y0, 'output_data\*.vtk'));

            perturb_contact = 0; direct_contact = 0;
            for ii = 1:length(fileinfo)
                snapshot = readVTK(fullfile(fileinfo(ii).folder, fileinfo(ii).name));
                XY = snapshot.points(:, 1:2);

                if min(pdist2(XY,[windward_edge_x, windward_edge_y],'euclidean','Smallest',1)) < fiber_beads_size  * 1.25 % if it's direct contact.

                    L_current = sqrt((XY(1,1)-XY(end,1))^2 + (XY(1,2)-XY(end,2))^2);
                    com_XY = mean(XY, 1); % length and CoM of current fiber

                    Gyr = 1/size(XY,1) * [sum((XY(:, 1)-com_XY(1)).^2),  sum((XY(:, 1)-com_XY(1)) .* (XY(:, 2)-com_XY(2)));
                        sum((XY(:, 2)-com_XY(2)) .* (XY(:, 1)-com_XY(1))), sum((XY(:, 2)-com_XY(2)).^2)];

                    [eigenV,eigenD] = eig(Gyr);
                    [d,ind] = sort(diag(eigenD));
                    Ds = eigenD(ind,ind);
                    Vs = eigenV(:,ind);
                    theta_c = atan(Vs(2,2)/Vs(1,2))/pi*180; % calculate theta_c
                    
                    L_start = [com_XY(1) - cosd(theta_c) * L_current, com_XY(2) - sind(theta_c)  * L_current];
                    L_end = [com_XY(1) + cosd(theta_c) * L_current, com_XY(2) + sind(theta_c)  * L_current];

                    P_inter = InterX([L_start_obsedge; L_end_obsedge]', [L_start; L_end]'); 
                    if ~isempty(P_inter)
                        y_c = (P_inter(2)- windward_edge_y(end)) / abs(windward_edge_y(1) - windward_edge_y(end)); % calculate y_c
                    else
                        y_c = nan;
                    end

                    break 
                end

            end

            if exist('theta_c','var')
                tmp_index =  thetheta0==deg_num & theL==L_num & they0==y0_num;
                Ind = find(tmp_index==1);
                Loc = ['I', num2str(Ind+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
                writematrix(y_c,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value inti the excel.
                Loc = ['J', num2str(Ind+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
                writematrix(theta_c,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value inti the excel.
            end

            clearvars theta_c y_c

        end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% plot: initial vs. contact %%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
xlsfile = readcell(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'results_2023_01_24.xlsx'],'Sheet','Sheet1','NumHeaderLines',1);
mask = cellfun(@ismissing, xlsfile); xlsfile(mask) = {nan};
together_plot = [cell2mat(xlsfile(:, 1:3)), cell2mat(xlsfile(:, 7)), ...
    cell2mat(xlsfile(:, 9:10)), cell2mat(xlsfile(:, 4))]; 
together_plot = together_plot';
prompt = {'The lower bound of the initial angle:', 'The upper bound of the initial angle:', ...
    'The lower bound of the contour length:','The upper bound of the contour length:'...
    'The lower bound of the initial position:','The upper bound of the initial position:'};
definput = {'-10', '10', 'nan', 'nan', 'nan', 'nan'};
answer = inputdlg(prompt, 'Input (please input NaN if there is no bound)', [1 35] , definput);

%%%%%%%%%%%%%%%%%%%%%%%%% assign the values %%%%%%%%%%%%%%%%%%%%%%%%%%
range_chi0_low = str2double(answer{1,1}); if isnan(range_chi0_low); range_chi0_low = -91; end
range_chi0_up = str2double(answer{2,1}); if isnan(range_chi0_up); range_chi0_up = 91; end
range_L_low = str2double(answer{3,1}); if isnan(range_L_low); range_L_low = 0; end
range_L_up = str2double(answer{4,1}); if isnan(range_L_up); range_L_up = 10; end
range_y0_low = str2double(answer{5,1}); if isnan(range_y0_low); range_y0_low = -10; end
range_y0_up = str2double(answer{6,1}); if isnan(range_y0_up); range_y0_up = 10; end
% chi_0
together_plot(:, together_plot(1, :) < range_chi0_low) = []; 
together_plot(:, together_plot(1, :) > range_chi0_up) = [];
% L
together_plot(:, together_plot(2, :) < range_L_low) = [];
together_plot(:, together_plot(2, :) > range_L_up) = [];
% y_0
together_plot(:, together_plot(3, :) < range_y0_low) = [];  
together_plot(:, together_plot(3, :) > range_y0_up) = [];

bypass_edge_together = together_plot(:, together_plot(4, :)==0); 
bypass_tip_together = together_plot(:, together_plot(4, :)==1);
pole_vaulting_together = together_plot(:, together_plot(4, :)==2); 
trapped_together = together_plot(:, together_plot(4, :)==3);



%%%%%%%%%%%%%%% plot theta_c vs y_c (with dynamics) %%%%%%%%%%%%%%%%%

figure('color', 'w'); set(gcf, 'Position', [100 100 1500 300]);
% for legend
plot(nan, nan, 'diamond','MarkerSize', 5,'MarkerEdgeColor','k','MarkerFaceColor','yellow'); 
plot(nan, nan, 'o','MarkerSize', 5,'MarkerEdgeColor','k','MarkerFaceColor','red');
plot(nan, nan,  '.','MarkerSize', 20,'MarkerEdgeColor', [0 .5 0]); 
plot(nan, nan,  '^','MarkerSize', 5,'MarkerEdgeColor','k','MarkerFaceColor','blue');  
% for plot
plot(trapped_together(6, :), trapped_together(5, :), 'diamond','MarkerSize', 5,'MarkerEdgeColor','k','MarkerFaceColor','yellow'); hold on 
plot(bypass_edge_together(6, :), bypass_edge_together(5, :), 'o','MarkerSize', 5,'MarkerEdgeColor','k','MarkerFaceColor','red'); hold on
plot(bypass_tip_together(6, :), bypass_tip_together(5, :),  '.','MarkerSize', 20,'MarkerEdgeColor', [0 .5 0]); hold on
plot(pole_vaulting_together(6, :), pole_vaulting_together(5, :),  '^','MarkerSize', 5,'MarkerEdgeColor','k','MarkerFaceColor','blue'); hold on

xlabel('$\theta_c$','FontSize', 18,'Interpreter', 'latex'); 
ylabel('$y_c$','FontSize', 18,'Interpreter', 'latex');
title_txt = ['$-10 < \theta_0 < 10$'];
title(title_txt,'FontSize', 18,'Interpreter', 'latex');
xlim([-90 90]); ylim([-0.1 1.1]);
legend({'Trapping','Below','Above','Pole-vaulting'}, 'Location', 'northwest','FontSize', 14,'Interpreter', 'latex')

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\theta_0m10to10_theta_c-y_c_dyn.png'],'Resolution',100)



%%%%%%%%%%%% plot theta_c vs y_c (with initial position y_0) %%%%%%%%%%%%%%

figure('color', 'w'); set(gcf, 'Position', [100 100 1500 300]);
% for legend
scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond'); hold on  % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'square'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', '^'); hold on % for legend only

scatter(trapped_together(6, :), trapped_together(5, :), 150, trapped_together(3, :), 'Filled', 'diamond','MarkerEdgeColor','k'); hold on 
scatter(bypass_edge_together(6, :), bypass_edge_together(5, :), 130, bypass_edge_together(3, :), 'Filled','MarkerEdgeColor','k'); hold on
scatter(bypass_tip_together(6, :), bypass_tip_together(5, :), 130, bypass_tip_together(3, :), 'Filled', 'square','MarkerEdgeColor','k'); hold on
scatter(pole_vaulting_together(6, :), pole_vaulting_together(5, :), 150, pole_vaulting_together(3, :), 'Filled', '^','MarkerEdgeColor','k'); hold on

% cmap(size(together_plot,2)); 
hcb=colorbar; caxis([0 1]); colormap jet
title(hcb,'$Initial\ position\ (y_0/h_{obs})$','FontSize', 16,'Interpreter', 'latex'); grid on

xlabel('$\theta_c$','FontSize', 18,'Interpreter', 'latex'); 
ylabel('$y_c$','FontSize', 18,'Interpreter', 'latex');
title_txt = ['$-10 < \theta_0 < 10$'];
title(title_txt,'FontSize', 18,'Interpreter', 'latex');
xlim([-90 90]); ylim([-0.1 1.1]);
legend({'Trapping','Below','Above','Pole-vaulting'}, 'Location', 'northwest','FontSize', 14,'Interpreter', 'latex')

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\theta_0m10to10_theta_c-y_c_y_0.png'],'Resolution',100)



%%%%%%%%%%%% plot theta_c vs y_c (with initial position theta_0) %%%%%%%%%%%%%%

figure('color', 'w'); set(gcf, 'Position', [100 100 1500 300]);
% for legend
scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond'); hold on  % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'square'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', '^'); hold on % for legend only

scatter(trapped_together(6, :), trapped_together(5, :), 150, trapped_together(1, :), 'Filled', 'diamond','MarkerEdgeColor','k'); hold on 
scatter(bypass_edge_together(6, :), bypass_edge_together(5, :), 130, bypass_edge_together(1, :), 'Filled','MarkerEdgeColor','k'); hold on
scatter(bypass_tip_together(6, :), bypass_tip_together(5, :), 130, bypass_tip_together(1, :), 'Filled', 'square','MarkerEdgeColor','k'); hold on
scatter(pole_vaulting_together(6, :), pole_vaulting_together(5, :), 150, pole_vaulting_together(1, :), 'Filled', '^','MarkerEdgeColor','k'); hold on

% cmap(size(together_plot,2)); 
hcb=colorbar; caxis([-10 10]); colormap jet
title(hcb,'$Initial\ angle\ ({\theta}_0)$','FontSize', 16,'Interpreter', 'latex'); grid on

xlabel('$\theta_c$','FontSize', 18,'Interpreter', 'latex'); 
ylabel('$y_c$','FontSize', 18,'Interpreter', 'latex');
title_txt = ['$-10 < \theta_0 < 10$'];
title(title_txt,'FontSize', 18,'Interpreter', 'latex');
xlim([-90 90]); ylim([-0.1 1.1]);
legend({'Trapping','Below','Above','Pole-vaulting'}, 'Location', 'northwest','FontSize', 14,'Interpreter', 'latex')

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\theta_0m10to10_theta_c-y_c_theta_0.png'],'Resolution',100)



%%%%%%%%%%%%%%% plot theta_c vs y_c (with contour length L) %%%%%%%%%%%%%%%%%

figure('color', 'w'); set(gcf, 'Position', [100 100 1500 300]);
% for legend
scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond'); hold on  % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'square'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', '^'); hold on % for legend only

scatter(trapped_together(6, :), trapped_together(5, :), 150, trapped_together(2, :), 'Filled', 'diamond','MarkerEdgeColor','k'); hold on 
scatter(bypass_edge_together(6, :), bypass_edge_together(5, :), 130, bypass_edge_together(2, :), 'Filled','MarkerEdgeColor','k'); hold on
scatter(bypass_tip_together(6, :), bypass_tip_together(5, :), 130, bypass_tip_together(2, :), 'Filled', 'square','MarkerEdgeColor','k'); hold on
scatter(pole_vaulting_together(6, :), pole_vaulting_together(5, :), 150, pole_vaulting_together(2, :), 'Filled', '^','MarkerEdgeColor','k'); hold on

% cmap(size(together_plot,2)); 
hcb=colorbar; caxis([0.5 1.5]); colormap jet
title(hcb,'$Contour\ length\ (L/l_{obs})$','FontSize', 16,'Interpreter', 'latex'); grid on

xlabel('$\theta_c$','FontSize', 18,'Interpreter', 'latex'); 
ylabel('$y_c$','FontSize', 18,'Interpreter', 'latex');
title_txt = ['$-10 < \theta_0 < 10$'];
title(title_txt,'FontSize', 18,'Interpreter', 'latex');
xlim([-90 90]); ylim([-0.1 1.1]);
legend({'Trapping','Below','Above','Pole-vaulting'}, 'Location', 'northwest','FontSize', 14,'Interpreter', 'latex')

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\theta_0m10to10_theta_c-y_c_L.png'],'Resolution',100)



%%%%%%%%%%%%%%% plot theta_c vs y_c (with deviation delta) %%%%%%%%%%%%%%%%%

figure('color', 'w'); set(gcf, 'Position', [100 100 1500 300]);
% for legend
scatter(nan, nan, 1, nan, 'filled', 'k', 'diamond'); hold on  % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', 'square'); hold on % for legend only
scatter(nan, nan, 1, nan, 'filled', 'k', '^'); hold on % for legend only

scatter(trapped_together(6, :), trapped_together(5, :), 150, trapped_together(7, :), 'Filled', 'diamond','MarkerEdgeColor','k'); hold on 
scatter(bypass_edge_together(6, :), bypass_edge_together(5, :), 130, bypass_edge_together(7, :), 'Filled','MarkerEdgeColor','k'); hold on
scatter(bypass_tip_together(6, :), bypass_tip_together(5, :), 130, bypass_tip_together(7, :), 'Filled', 'square','MarkerEdgeColor','k'); hold on
scatter(pole_vaulting_together(6, :), pole_vaulting_together(5, :), 150, pole_vaulting_together(7, :), 'Filled', '^','MarkerEdgeColor','k'); hold on

% cmap(size(together_plot,2)); 
hcb=colorbar; caxis([-0.6 0.6]); colormap jet
title(hcb,'$Deviation\ (\delta)$','FontSize', 16,'Interpreter', 'latex'); grid on

xlabel('$\theta_c$','FontSize', 18,'Interpreter', 'latex'); 
ylabel('$y_c$','FontSize', 18,'Interpreter', 'latex');
title_txt = ['$-10 < \theta_0 < 10$'];
title(title_txt,'FontSize', 18,'Interpreter', 'latex');
xlim([-90 90]); ylim([-0.1 1.1]);
legend({'Trapping','Below','Above','Pole-vaulting'}, 'Location', 'northwest','FontSize', 14,'Interpreter', 'latex')

% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\theta_0m10to10_theta_c-y_c_delta.png'],'Resolution',100)



%%%%%%%%%%%%%%% plot (theta_0, y_0) vs (theta_c, y_c) in vector field %%%%%%%%%%%%%%%%%

% let the data 'align' for reshape
together_plot(:, together_plot(3, :) == 0.1875) = [];
together_plot(:, together_plot(3, :) == 0.3125) = [];
together_plot(:, together_plot(3, :) == 0.4375) = [];
together_plot(:, together_plot(3, :) == 0.21875) = [];
together_plot(:, together_plot(3, :) == 0.28125) = [];
together_plot(:, together_plot(3, :) == 0.34375) = [];
together_plot(:, together_plot(3, :) == 0.40625) = [];

% choose an 'L' to plot:
choose_L = 0.5;
together_plot(:, together_plot(2, :) ~= choose_L) = [];

toPlot_theta_0 = together_plot(1, :);
toPlot_L_0 = together_plot(2, :);
toPlot_y_0 = together_plot(3, :);
toPlot_y_c = together_plot(5, :);
toPlot_theta_c = together_plot(6, :);

[toPlot_y_0_X,toPlot_theta_0_Y] = meshgrid(unique(toPlot_theta_0),unique(toPlot_y_0));

toPlot_theta_c = reshape(toPlot_theta_c, [length(unique(toPlot_y_0)), length(unique(toPlot_theta_0))]);
toPlot_y_c = reshape(toPlot_y_c, [length(unique(toPlot_y_0)), length(unique(toPlot_theta_0))]);
toPlot_contact_X = toPlot_y_c .* sind(toPlot_theta_c);
toPlot_contact_Y = toPlot_y_c .* cosd(toPlot_theta_c);

toPlot_contact_X_base = ~isnan(toPlot_y_c) .* ones(size(toPlot_y_c,1), size(toPlot_y_c,2)) .* sind(toPlot_theta_c);
toPlot_contact_Y_base = ~isnan(toPlot_y_c) .* ones(size(toPlot_y_c,1), size(toPlot_y_c,2)) .* cosd(toPlot_theta_c);

figure('color', 'w'); set(gcf, 'Position', [100 100 1500 300]);
quiver(toPlot_y_0_X,toPlot_theta_0_Y,toPlot_contact_Y_base,toPlot_contact_X_base, ...
    'k','LineWidth', 1, 'ShowArrowHead','off'); hold on
quiver(toPlot_y_0_X,toPlot_theta_0_Y,toPlot_contact_Y,toPlot_contact_X,'m','LineWidth', 2); axis equal
xlabel('$Initial\ angle\ (\theta_0)$','FontSize', 18,'Interpreter', 'latex'); 
ylabel('$Initial\ position\ (y_0)$','FontSize', 18,'Interpreter', 'latex');
title_txt = ['$Length\ L = ',num2str(choose_L), '$'];
title(title_txt,'FontSize', 18,'Interpreter', 'latex');
xlim([-12 12]); ylim([-0.1 2.1]); grid on
% f=gcf;
% exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
%     '\Figures\about contact information vs initial condition\contourL_', num2str(choose_L) ,'.png'],'Resolution',100)