clear; close all; clc;

[excelname, excelpathname] = uigetfile(['D:\Dropbox\Collaboration - LadHyX\' ...
    'Give_to_Zhibo_nonShared\.xlsx'], 'Please choose the excel accordingly and carefully!!');
xlsfile = readcell([excelpathname, excelname],'Sheet','Sheet1','NumHeaderLines',1); 
thedeg = cell2mat(xlsfile(:, 1)); 
theL = round(cell2mat(xlsfile(:, 2)), 1); 
they0 = round(cell2mat(xlsfile(:, 3)), 3); 

% obs = readVTK('D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\obstacle_beads.vtk');
% obs_beads_size = 1.5e-6; % notice here (radius)
fiber_beads_size = 0.8e-6 * 1.25; % notice here (radius)
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

                if_direct_contect = min(pdist2(XY,obs_2d,'euclidean','Smallest',1)) < fiber_beads_size; % if it's direct contact.
                if if_direct_contect
                    direct_contact = direct_contact + 1;
                end

%                 in_triangle = inpolygon(XY(:,1),XY(:,2),obs_2d(:,1),obs_2d(:,2)); % the points on the fiber which are inside the triangular obstacle
%                 in_perturb_C1 = inpolygon(XY(:,1),XY(:,2),pert_C1(:,1),pert_C1(:,2)); % the points on the fiber which are inside the perturbed line 1
%                 if_inbetween_twolines = and(~logical(in_triangle), logical(in_perturb_C1)); % the fiber is in between the two contours
%                 if_inbetween_twolines = logical(sum(if_inbetween_twolines));
% 
%                 if min(pdist2(XY,pert_C1,'euclidean','Smallest',1)) < fiber_beads_size || ...
%                     if_direct_contect || if_inbetween_twolines
%                     perturb_contact = perturb_contact + 1;
%                 end

                in_perturb_C2 = inpolygon(XY(:,1),XY(:,2),pert_C2(:,1),pert_C2(:,2)); % the points on the fiber which are inside the perturbed line 2
                if_in_C2 = logical(sum(logical(in_perturb_C2)));

                if (min(pdist2(XY,pert_C2,'euclidean','Smallest',1)) < fiber_beads_size || if_in_C2) && ~if_direct_contect % counting number of perturbed contact
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
