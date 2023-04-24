%%% To calculate the evolution of the fiber orientation
% starts from zero degree and varies continuously.
% might be useful in the future

clear; clc;
commandwindow;

load(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\input_data\obs_2d_20230223.mat']); 
obs_x = mean(obs_2d(:, 1)); h_obs = max(obs_2d(:, 2)) - min(obs_2d(:, 2));

parent_path = 'D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\Data_Give_to_Zhibo_20230223\simulations';
sub1_path = dir(parent_path);
for sub1Path_i = 3:length(sub1_path)
    current_deg = sub1_path(sub1Path_i).name; 
    newStr = strrep(current_deg,'o','.'); newStr = strrep(newStr,'m','-'); 
    deg_num = str2double(newStr(8:end));
    sub2_path = dir(fullfile(parent_path, current_deg));
    for sub2Path_i = 3:length(sub2_path)
        current_L = sub2_path(sub2Path_i).name;
        sub3_path = dir(fullfile(parent_path, current_deg, current_L));
        for sub3Path_i = 3:length(sub3_path)
            current_y0 = sub3_path(sub3Path_i).name;
            fileinfo = dir(fullfile(parent_path, current_deg, current_L, current_y0, 'output_data\*.vtk'));

            x = nan(1, length(fileinfo)); 
            snapshot = readVTK(fullfile(fileinfo(1).folder, fileinfo(1).name));
            x(1) = (mean(snapshot.points(:, 1)) - obs_x) / h_obs;

            ori_ee = nan(1, length(fileinfo));
            if deg_num >= 0
                ori_ee(1) = deg_num;
            else
                ori_ee(1) = 360 + deg_num;  % convert initial angle to range of [0 ,360)
            end
            for ii = 2:length(fileinfo)
                snapshot = readVTK(fullfile(fileinfo(ii).folder, fileinfo(ii).name));
                x(ii) = (mean(snapshot.points(:, 1)) - obs_x) / h_obs;

                XY_1 = snapshot.points(1, 1:2);
                XY_2 = snapshot.points(end, 1:2);
                %     d_XY = XY_2-XY_1;                         % Difference
                %     quiver(XY_1(1),XY_2(2),d_XY(1),d_XY(2),0); hold on

                if XY_2(1)-XY_1(1) > 0 && XY_2(2)-XY_1(2) >= 0  % the first quadrant
                    ori_ee(ii) = atand((XY_2(2)-XY_1(2)) / (XY_2(1)-XY_1(1)));  % ori_ee: [0, 360); X positive is 0°.
                elseif XY_2(1)-XY_1(1) <= 0 && XY_2(2)-XY_1(2) > 0  % the second quadrant
                    ori_ee(ii) = 180 + atand((XY_2(2)-XY_1(2)) / (XY_2(1)-XY_1(1)));
                elseif XY_2(1)-XY_1(1) < 0 && XY_2(2)-XY_1(2) <= 0  % the third quadrant
                    ori_ee(ii) = 180 + atand((XY_2(2)-XY_1(2)) / (XY_2(1)-XY_1(1)));
                elseif XY_2(1)-XY_1(1) >= 0 && XY_2(2)-XY_1(2) < 0  % the fourth quadrant
                    ori_ee(ii) = 360 + atand((XY_2(2)-XY_1(2)) / (XY_2(1)-XY_1(1)));
                end
            end

            %%%%% To make the curve continuous and start from zero
            ori_ee = ori_ee - ori_ee(1); 
            ori_ee = [0, ori_ee];
            [pks1, loc1] = findpeaks(diff(ori_ee), 'MinPeakHeight', 180);
            [pks2, loc2] = findpeaks(-diff(ori_ee), 'MinPeakHeight', 180);

            if ~isempty(loc1) || ~isempty(loc2)
                loc = sort([loc1, loc2]);
                if logical(sum(ismember(loc(1), loc1)))
                    for foo = 1:length(loc1)
                        try
                            ori_ee(loc1(foo)+1: loc2(foo)) = ori_ee(loc1(foo)+1: loc2(foo)) - 360;
                        catch
                            ori_ee(loc1(foo)+1: end) = ori_ee(loc1(foo)+1:end) - 360;
                        end
                    end
                else
                    for foo = 1:length(loc2)
                        try
                            ori_ee(loc2(foo)+1: loc1(foo)) = ori_ee(loc2(foo)+1: loc1(foo)) + 360;
                        catch
                            ori_ee(loc2(foo)+1: end) = ori_ee(loc2(foo)+1: end) + 360;
                        end
                    end
                end
            end
            ori_ee(1) = [];
%             ori_ee(ori_ee>180) = ori_ee(ori_ee>180) - 360;
            
            figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
            plot(ori_ee, 'o','MarkerSize', 8,'MarkerEdgeColor','k','MarkerFaceColor','red');

            xlabel('$Frame$','FontSize', 18,'Interpreter', 'latex');
            ylabel('$\Delta\theta$','FontSize', 18,'Interpreter', 'latex');

            set_plot(gcf, gca)
            f=gcf;
            exportgraphics(f,['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared' ...
                '\Data_Give_to_Zhibo_20230223_videos\Orientations\', current_deg, '_', current_L, ...
                '_', current_y0, '_animation.png'],'Resolution', 100)
            close all

            save(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared' ...
                '\Data_Give_to_Zhibo_20230223_videos\Orientations\', current_deg, '_', current_L, ...
                '_', current_y0, '_orientations.mat'], 'ori_ee', 'x');
        end
    end
end


%% orientations vs x: given theta_0 and L to plot different y_0
clear; close all; clc;

Files = dir(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared' ...
    '\Data_Give_to_Zhibo_20230223_videos\Orientations\*.mat']);

for choose_theta0 = -10:2.5:10
    for choose_L = [0.5:0.1:1, 1.2:0.2:1.4]

        n = 1;
        for ii = 1: length(Files)

            load(fullfile(Files(ii).folder, Files(ii).name));
            current_deg =  extractBetween(Files(ii).name,'theta0_', '_L');
            current_deg = strrep(current_deg,'o','.'); current_deg = strrep(current_deg,'m','-');
            current_L =  extractBetween(Files(ii).name,'L_', '_y0');
            current_L = strrep(current_L,'o','.');
            current_y0 =  extractBetween(Files(ii).name,'y0_', '_ori');
            current_y0 = strrep(current_y0,'o','.');
            if isnan(str2double(current_y0{1}))
                current_y0_tmp = current_y0{1};
                current_y0{1} = current_y0_tmp(1:end-10);
            end

            if str2double(current_deg) ~= choose_theta0 || str2double(current_L) ~= choose_L
                % choose the y_0 and L to be plotted
                continue
            else
                title_deg = current_deg;
                title_L = current_L;
                To_Plot{n, 1} = ori_ee;
                To_Plot{n, 2} = x;
                To_Plot{n, 3} = str2double(current_y0); % to be sorted and plotted
                n = n + 1;
            end

        end

        figure('color', 'w'); set(gcf, 'Position', [100 100 1000 500]);
        cmap = cmocean('thermal');  legend_txt = {};
        [~, sortID] = sort(cell2mat(To_Plot(:,3)));  % sort the plotting order

        color_ind = 1;
        for jj = [3 5 6 7 8 9 13 16]%1:size(To_Plot, 1)

            plot(To_Plot{sortID(jj), 2}, To_Plot{sortID(jj), 1}, 'o','MarkerSize', 5,'MarkerEdgeColor','k', ...
                'MarkerFaceColor', cmap(color_ind*30,:)); hold on

            xlabel('$x/h_{obs}$','FontSize', 18,'Interpreter', 'latex');
            ylabel('$\Delta\theta$','FontSize', 18,'Interpreter', 'latex');

            title(strcat('$\theta_0=', title_deg, '^{\circ}\ and\ L=', title_L, '$'), ...
                'FontSize', 24, 'Interpreter', 'latex');
            legend_txt = [legend_txt, strcat('$y_0=', num2str(To_Plot{sortID(jj), 3}),'$')];

            color_ind = color_ind + 1;

        end
        xlim([-10 10])
        ax = gca; ax.FontSize = 18;
        legend(legend_txt, 'Location', 'northwest', 'FontSize', 18,  'Interpreter', 'latex');

        savename = strcat('theta_0=', title_deg, '_L=', title_L, '_orientation_vs_x.png');

        set_plot(gcf, gca)
        f=gcf;
        exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\' ...
            'Figures\about orientations vs x\Given theta0 and L\', savename{1}],'Resolution',100)

        close
    end
end



%% orientations vs x: given y_0 and L to plot different theta_0
clear; close all; clc;

Files = dir(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared' ...
    '\Data_Give_to_Zhibo_20230223_videos\Orientations\*.mat']);

for choose_y0 = [0.1:0.1:0.9, 0.15:0.1:0.55, 0.325:0.05:0.525]
    for choose_L = [0.5:0.1:1, 1.2:0.2:1.4]

        n = 1;
        for ii = 1: length(Files)

            load(fullfile(Files(ii).folder, Files(ii).name));
            current_deg =  extractBetween(Files(ii).name,'theta0_', '_L');
            current_deg = strrep(current_deg,'o','.'); current_deg = strrep(current_deg,'m','-');
            current_L =  extractBetween(Files(ii).name,'L_', '_y0');
            current_L = strrep(current_L,'o','.');
            current_y0 =  extractBetween(Files(ii).name,'y0_', '_ori');
            current_y0 = strrep(current_y0,'o','.');
            if isnan(str2double(current_y0{1}))
                current_y0_tmp = current_y0{1};
                current_y0{1} = current_y0_tmp(1:end-10);
            end

            if str2double(current_y0) ~= choose_y0 || str2double(current_L) ~= choose_L
                % choose the y_0 and L to be plotted
                continue
            else
                title_y0 = current_y0;
                title_L = current_L;
                To_Plot{n, 1} = ori_ee;
                To_Plot{n, 2} = x;
                To_Plot{n, 3} = str2double(current_deg); % to be sorted and plotted
                n = n + 1;
            end

        end

        figure('color', 'w'); set(gcf, 'Position', [100 100 1000 500]);
        cmap = cmocean('thermal');  legend_txt = {};
        [~, sortID] = sort(cell2mat(To_Plot(:,3)));  % sort the plotting order

        color_ind = 1;
        for jj = 1:size(To_Plot, 1)

            plot(To_Plot{sortID(jj), 2}, To_Plot{sortID(jj), 1}, 'o','MarkerSize', 5,'MarkerEdgeColor','k', ...
                'MarkerFaceColor', cmap(color_ind*floor(255/size(To_Plot, 1)),:)); hold on

            xlabel('$x/h_{obs}$','FontSize', 18,'Interpreter', 'latex');
            ylabel('$\Delta\theta$','FontSize', 18,'Interpreter', 'latex');

            title(strcat('$y_0=', title_y0, '\ and\ L=', title_L, '$'), ...
                'FontSize', 24, 'Interpreter', 'latex');
            legend_txt = [legend_txt, strcat('$\theta_0=', num2str(To_Plot{sortID(jj), 3}),'^{\circ}$')];

            color_ind = color_ind + 1;

        end
        xlim([-10 10])
        ax = gca; ax.FontSize = 18;
        legend(legend_txt, 'Location', 'northwest', 'FontSize', 18,  'Interpreter', 'latex');

        savename = strcat('y_0=', title_y0, '_L=', title_L, '_orientation_vs_x.png');

        set_plot(gcf, gca)
        f=gcf;
        exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\' ...
            'Figures\about orientations vs x\Given y0 and L\', savename{1}],'Resolution',100)

        close
    end
end



%% orientations vs x: given y_0 and theta_0 to plot different L
clear; close all; clc;

Files = dir(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared' ...
    '\Data_Give_to_Zhibo_20230223_videos\Orientations\*.mat']);

for choose_y0 = [0.1:0.1:0.9, 0.15:0.1:0.55, 0.325:0.05:0.525]
    for choose_theta0 = -10:2.5:10

        n = 1;
        for ii = 1: length(Files)

            load(fullfile(Files(ii).folder, Files(ii).name));
            current_deg =  extractBetween(Files(ii).name,'theta0_', '_L');
            current_deg = strrep(current_deg,'o','.'); current_deg = strrep(current_deg,'m','-');
            current_L =  extractBetween(Files(ii).name,'L_', '_y0');
            current_L = strrep(current_L,'o','.');
            current_y0 =  extractBetween(Files(ii).name,'y0_', '_ori');
            current_y0 = strrep(current_y0,'o','.');
            if isnan(str2double(current_y0{1}))
                current_y0_tmp = current_y0{1};
                current_y0{1} = current_y0_tmp(1:end-10);
            end

            if str2double(current_y0) ~= choose_y0 || str2double(current_deg) ~= choose_theta0
                % choose the y_0 and L to be plotted
                continue
            else
                title_y0 = current_y0;
                title_deg = current_deg;
                To_Plot{n, 1} = ori_ee;
                To_Plot{n, 2} = x;
                To_Plot{n, 3} = str2double(current_L); % to be sorted and plotted
                n = n + 1;
            end

        end

        figure('color', 'w'); set(gcf, 'Position', [100 100 1000 500]);
        cmap = cmocean('thermal');  legend_txt = {};
        [~, sortID] = sort(cell2mat(To_Plot(:,3)));  % sort the plotting order

        color_ind = 1;
        for jj = 1:size(To_Plot, 1)

            plot(To_Plot{sortID(jj), 2}, To_Plot{sortID(jj), 1}, 'o','MarkerSize', 5,'MarkerEdgeColor','k', ...
                'MarkerFaceColor', cmap(color_ind*floor(255/size(To_Plot, 1)),:)); hold on

            xlabel('$x/h_{obs}$','FontSize', 18,'Interpreter', 'latex');
            ylabel('$\Delta\theta$','FontSize', 18,'Interpreter', 'latex');

            title(strcat('$y_0=', title_y0, '\ and\ \theta_0=', title_deg, '^{\circ}$'), ...
                'FontSize', 24, 'Interpreter', 'latex');
            legend_txt = [legend_txt, strcat('$L=', num2str(To_Plot{sortID(jj), 3}),'$')];

            color_ind = color_ind + 1;

        end
        xlim([-10 10])
        ax = gca; ax.FontSize = 18;
        legend(legend_txt, 'Location', 'northwest', 'FontSize', 18,  'Interpreter', 'latex');

        savename = strcat('y_0=', title_y0, '_theta_0=', title_deg, '_orientation_vs_x.png');

        set_plot(gcf, gca)
        f=gcf;
        exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\' ...
            'Figures\about orientations vs x\Given y0 and theta0\', savename{1}],'Resolution',100)

        close
    end
end



%% orientations vs x: plot until contact
clear; close all; clc;

xlsfile = readcell(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'Data_Give_to_Zhibo_20230223\results_2023_02_23_Zhibo.xlsx'], ...
    'Sheet','Sheet1','NumHeaderLines',1);
mask = cellfun(@ismissing, xlsfile); xlsfile(mask) = {nan};
thedeg = cell2mat(xlsfile(:, 1)); 
theL = round(cell2mat(xlsfile(:, 2)), 1); 
they0 = round(cell2mat(xlsfile(:, 3)), 3); 

ite_contact = cell2mat(xlsfile(:, 11)); 


Files = dir(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared' ...
    '\Data_Give_to_Zhibo_20230223_videos\Orientations\*.mat']);

for ii = 1: length(Files)

    load(fullfile(Files(ii).folder, Files(ii).name));

    savename = strcat(Files(ii).name(1:end-4), '_untilContact.png');

    current_deg =  extractBetween(Files(ii).name,'theta0_', '_L');
    newStr = strrep(current_deg,'o','.'); newStr = strrep(newStr,'m','-');
    deg_num = str2double(newStr);
    current_L =  extractBetween(Files(ii).name,'L_', '_y0');
    newStr = strrep(current_L,'o','.');
    L_num = str2double(newStr);
    current_y0 =  extractBetween(Files(ii).name,'y0_', '_ori');
    newStr = strrep(current_y0,'o','.');
    if isnan(str2double(newStr{1}))
        current_y0_tmp = newStr{1};
        newStr{1} = current_y0_tmp(1:end-10);
        y0_num = str2double(newStr);
    else
        y0_num = str2double(newStr);
    end

    tmp_index =  thedeg==deg_num & theL==L_num & ismembertol(they0, y0_num,0.0125);
    Ind = find(tmp_index==1);

    if ~isnan(ite_contact(Ind))
        figure('color', 'w'); set(gcf, 'Position', [100 100 1000 500]);
        x_untilContact = x(36:ite_contact(Ind));
        ori_ee_untilContact = ori_ee(36:ite_contact(Ind));
        plot(x_untilContact, ori_ee_untilContact, 'o','MarkerSize', 10,'MarkerEdgeColor','k', ...
            'MarkerFaceColor', 'm'); 

        xlabel('$x/h_{obs}$','FontSize', 18,'Interpreter', 'latex');
        ylabel('$\Delta\theta$','FontSize', 18,'Interpreter', 'latex');
        xlim([-5 2])
        ax = gca; ax.FontSize = 18;

        set_plot(gcf, gca)
        f=gcf;
        exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle' ...
            '\Figures\about orientations vs x\until contact\', savename],'Resolution',100)

        save(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
            'Data_Give_to_Zhibo_20230223_videos\Orientations\until contact\theta0_', ...
            current_deg{1}, '_L_', current_L{1},'_y0_', current_y0{1}, '_orientations_untilContact.mat'], ...
            'ori_ee', 'x', 'x_untilContact', 'ori_ee_untilContact');

    end

    close
end



%% orientations vs x: given theta_0 and L to plot different y_0 (until contact)
clear; close all; clc;

Files = dir(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared' ...
    '\Data_Give_to_Zhibo_20230223_videos\Orientations\until contact\*.mat']);

for choose_theta0 = -10:2.5:10
    for choose_L = [0.5:0.1:1, 1.2:0.2:1.4]

        n = 1;
        for ii = 1: length(Files)

            load(fullfile(Files(ii).folder, Files(ii).name));
            current_deg =  extractBetween(Files(ii).name,'theta0_', '_L');
            current_deg = strrep(current_deg,'o','.'); current_deg = strrep(current_deg,'m','-');
            current_L =  extractBetween(Files(ii).name,'L_', '_y0');
            current_L = strrep(current_L,'o','.');
            current_y0 =  extractBetween(Files(ii).name,'y0_', '_ori');
            current_y0 = strrep(current_y0,'o','.');
            if isnan(str2double(current_y0{1}))
                current_y0_tmp = current_y0{1};
                current_y0{1} = current_y0_tmp(1:end-10);
            end

            if str2double(current_deg) ~= choose_theta0 || str2double(current_L) ~= choose_L
                % choose the y_0 and L to be plotted
                continue
            else
                title_deg = current_deg;
                title_L = current_L;
                To_Plot{n, 1} = ori_ee_untilContact;
                To_Plot{n, 2} = x_untilContact;
                To_Plot{n, 3} = str2double(current_y0); % to be sorted and plotted
                n = n + 1;
            end

        end

        figure('color', 'w'); set(gcf, 'Position', [100 100 1000 500]);
        cmap = cmocean('thermal');  legend_txt = {};
        [~, sortID] = sort(cell2mat(To_Plot(:,3)));  % sort the plotting order

        color_ind = 1;
        for jj = 1:size(To_Plot, 1)

            plot(To_Plot{sortID(jj), 2}, To_Plot{sortID(jj), 1}, 'o','MarkerSize', 5,'MarkerEdgeColor','k', ...
                'MarkerFaceColor', cmap(color_ind*floor(255/size(To_Plot, 1)),:)); hold on

            xlabel('$x/h_{obs}$','FontSize', 18,'Interpreter', 'latex');
            ylabel('$\Delta\theta$','FontSize', 18,'Interpreter', 'latex');

            title(strcat('$\theta_0=', title_deg, '^{\circ}\ and\ L=', title_L, '$'), ...
                'FontSize', 24, 'Interpreter', 'latex');
            legend_txt = [legend_txt, strcat('$y_0=', num2str(To_Plot{sortID(jj), 3}),'$')];

            color_ind = color_ind + 1;

        end
        xlim([-5 0])
        ax = gca; ax.FontSize = 18;
        legend(legend_txt, 'Location', 'northwest', 'FontSize', 18,  'Interpreter', 'latex');

        savename = strcat('theta_0=', title_deg, '_L=', title_L, '_orientation_vs_x.png');

        set_plot(gcf, gca)
        f=gcf;
        exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\' ...
            'Figures\about orientations vs x\until contact\Given theta0 and L\', savename{1}],'Resolution',100)

        close
    end
end



%% orientations vs x: given y_0 and L to plot different theta_0 (until contact)
clear; close all; clc;

Files = dir(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared' ...
    '\Data_Give_to_Zhibo_20230223_videos\Orientations\until contact\*.mat']);

for choose_y0 = [0.1:0.1:0.9, 0.15:0.1:0.55, 0.325:0.05:0.525]
    for choose_L = [0.5:0.1:1, 1.2:0.2:1.4]

        n = 1;
        for ii = 1: length(Files)

            if ~exist(fullfile(Files(ii).folder, Files(ii).name), 'file')
                continue
            else
                load(fullfile(Files(ii).folder, Files(ii).name));
                current_deg =  extractBetween(Files(ii).name,'theta0_', '_L');
                current_deg = strrep(current_deg,'o','.'); current_deg = strrep(current_deg,'m','-');
                current_L =  extractBetween(Files(ii).name,'L_', '_y0');
                current_L = strrep(current_L,'o','.');
                current_y0 =  extractBetween(Files(ii).name,'y0_', '_ori');
                current_y0 = strrep(current_y0,'o','.');
                if isnan(str2double(current_y0{1}))
                    current_y0_tmp = current_y0{1};
                    current_y0{1} = current_y0_tmp(1:end-10);
                end

                if str2double(current_y0) ~= choose_y0 || str2double(current_L) ~= choose_L
                    % choose the y_0 and L to be plotted
                    continue
                else
                    title_y0 = current_y0;
                    title_L = current_L;
                    To_Plot{n, 1} = ori_ee_untilContact;
                    To_Plot{n, 2} = x_untilContact;
                    To_Plot{n, 3} = str2double(current_deg); % to be sorted and plotted
                    n = n + 1;
                end

            end
        end
        
        figure('color', 'w'); set(gcf, 'Position', [100 100 1000 500]);
        cmap = cmocean('thermal');  legend_txt = {};
        [~, sortID] = sort(cell2mat(To_Plot(:,3)));  % sort the plotting order

        color_ind = 1;
        for jj = 1:size(To_Plot, 1)

            plot(To_Plot{sortID(jj), 2}, To_Plot{sortID(jj), 1}, 'o','MarkerSize', 5,'MarkerEdgeColor','k', ...
                'MarkerFaceColor', cmap(color_ind*floor(255/size(To_Plot, 1)),:)); hold on

            xlabel('$x/h_{obs}$','FontSize', 18,'Interpreter', 'latex');
            ylabel('$\Delta\theta$','FontSize', 18,'Interpreter', 'latex');

            title(strcat('$y_0=', title_y0, '\ and\ L=', title_L, '$'), ...
                'FontSize', 24, 'Interpreter', 'latex');
            legend_txt = [legend_txt, strcat('$\theta_0=', num2str(To_Plot{sortID(jj), 3}),'^{\circ}$')];

            color_ind = color_ind + 1;

        end
        xlim([-5 0])
        ax = gca; ax.FontSize = 18;
        legend(legend_txt, 'Location', 'northwest', 'FontSize', 18,  'Interpreter', 'latex');

        savename = strcat('y_0=', title_y0, '_L=', title_L, '_orientation_vs_x.png');

        set_plot(gcf, gca)
        f=gcf;
        exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\' ...
            'Figures\about orientations vs x\until contact\Given y0 and L\', savename{1}],'Resolution',100)

        close
    end
end



%% orientations vs x: given y_0 and theta_0 to plot different L (until contact)
clear; close all; clc;

Files = dir(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared' ...
    '\Data_Give_to_Zhibo_20230223_videos\Orientations\until contact\*.mat']);

for choose_y0 = [0.1:0.1:0.9, 0.15:0.1:0.55, 0.325:0.05:0.525]
    for choose_theta0 = -10:2.5:10

        n = 1;
        for ii = 1: length(Files)

            if ~exist(fullfile(Files(ii).folder, Files(ii).name), 'file')
                continue
            else
                load(fullfile(Files(ii).folder, Files(ii).name));
                current_deg =  extractBetween(Files(ii).name,'theta0_', '_L');
                current_deg = strrep(current_deg,'o','.'); current_deg = strrep(current_deg,'m','-');
                current_L =  extractBetween(Files(ii).name,'L_', '_y0');
                current_L = strrep(current_L,'o','.');
                current_y0 =  extractBetween(Files(ii).name,'y0_', '_ori');
                current_y0 = strrep(current_y0,'o','.');
                if isnan(str2double(current_y0{1}))
                    current_y0_tmp = current_y0{1};
                    current_y0{1} = current_y0_tmp(1:end-10);
                end

                if str2double(current_y0) ~= choose_y0 || str2double(current_deg) ~= choose_theta0
                    % choose the y_0 and L to be plotted
                    continue
                else
                    title_y0 = current_y0;
                    title_deg = current_deg;
                    To_Plot{n, 1} = ori_ee_untilContact;
                    To_Plot{n, 2} = x_untilContact;
                    To_Plot{n, 3} = str2double(current_L); % to be sorted and plotted
                    n = n + 1;
                end

            end
        end

        figure('color', 'w'); set(gcf, 'Position', [100 100 1000 500]);
        cmap = cmocean('thermal');  legend_txt = {};
        To_Plot([1,3,5,7], :) = []; % only plot L=0.6, 0.8, 1 and 1.4.
        [~, sortID] = sort(cell2mat(To_Plot(:,3)));  % sort the plotting order

        color_ind = 1;
        for jj = 1:size(To_Plot, 1)

            plot(To_Plot{sortID(jj), 2}, To_Plot{sortID(jj), 1}, 'o','MarkerSize', 8,'MarkerEdgeColor','k', ...
                'MarkerFaceColor', cmap(color_ind*floor(255/size(To_Plot, 1)),:)); hold on

            xlabel('$x/h_{obs}$','FontSize', 18,'Interpreter', 'latex');
            ylabel('$\Delta\theta$','FontSize', 18,'Interpreter', 'latex');

            title(strcat('$y_0=', title_y0, '\ and\ \theta_0=', title_deg, '^{\circ}$'), ...
                'FontSize', 24, 'Interpreter', 'latex');
            legend_txt = [legend_txt, strcat('$L=', num2str(To_Plot{sortID(jj), 3}),'$')];

            color_ind = color_ind + 1;

        end
        xlim([-5 0])
        ax = gca; ax.FontSize = 18;
        legend(legend_txt, 'Location', 'northwest', 'FontSize', 18,  'Interpreter', 'latex');

        savename = strcat('y_0=', title_y0, '_theta_0=', title_deg, '_orientation_vs_x.png');

        set_plot(gcf, gca) 
        f=gcf;
        exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\' ...
            'Figures\about orientations vs x\until contact\Given y0 and theta0\', savename{1}],'Resolution',100)

        close
    end
end



%% orientations vs x: given y_0 (0.425) and theta_0 (-7.5) to plot different L (until contact)
clear; close all; clc;

Files = dir(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared' ...
    '\Data_Give_to_Zhibo_20230223_videos\Orientations\until contact\*.mat']);

for choose_y0 = 0.425
    for choose_theta0 = -7.5

        n = 1;
        for ii = 1: length(Files)

            if ~exist(fullfile(Files(ii).folder, Files(ii).name), 'file')
                continue
            else
                load(fullfile(Files(ii).folder, Files(ii).name));
                current_deg =  extractBetween(Files(ii).name,'theta0_', '_L');
                current_deg = strrep(current_deg,'o','.'); current_deg = strrep(current_deg,'m','-');
                current_L =  extractBetween(Files(ii).name,'L_', '_y0');
                current_L = strrep(current_L,'o','.');
                current_y0 =  extractBetween(Files(ii).name,'y0_', '_ori');
                current_y0 = strrep(current_y0,'o','.');
                if isnan(str2double(current_y0{1}))
                    current_y0_tmp = current_y0{1};
                    current_y0{1} = current_y0_tmp(1:end-10);
                end

                if str2double(current_y0) ~= choose_y0 || str2double(current_deg) ~= choose_theta0
                    % choose the y_0 and L to be plotted
                    continue
                else
                    title_y0 = current_y0;
                    title_deg = current_deg;
                    To_Plot{n, 1} = ori_ee_untilContact;
                    To_Plot{n, 2} = x_untilContact;
                    To_Plot{n, 3} = str2double(current_L); % to be sorted and plotted
                    n = n + 1;
                end

            end
        end

        fig1 = figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);
        cmap = cmocean('matter');  legend_txt = {};
        To_Plot([1,3,5,7], :) = []; % only plot L=0.6, 0.8, 1 and 1.4.
        [~, sortID] = sort(cell2mat(To_Plot(:,3)));  % sort the plotting order

        color_ind = size(To_Plot, 1);
        for jj = 1:size(To_Plot, 1)

            plot(To_Plot{sortID(jj), 2}, To_Plot{sortID(jj), 1}, 'o','MarkerSize', 8,'MarkerEdgeColor','k', ...
                'MarkerFaceColor', cmap(color_ind*60,:)); hold on

            xlabel('$x/h_{\rm obs}$','FontSize', 24,'Interpreter', 'latex');
            ylabel('$\Delta\theta$','FontSize', 24,'Interpreter', 'latex');

            color_ind = color_ind - 1; % coincidence with the vector plot

        end
        xlim([-5 0])
        set(gca, ...
            'Box', 'On', ...
            'XGrid', 'On', ...
            'YGrid', 'On', ...
            'GridAlpha', 0.5, ...
            'FontSize', 24, ...
            'NextPlot','replacechildren', ...
            'TickLabelInterpreter','latex')

        axes('Position',[.25 .3 .4 .3])

        color_ind = size(To_Plot, 1);
        for jj = 1:size(To_Plot, 1)

            select_To_Plot = To_Plot{sortID(jj), 2};
            select_ind = select_To_Plot < -4;
            plot(To_Plot{sortID(jj), 2}(select_ind), To_Plot{sortID(jj), 1}((select_ind)), ...
                'o','MarkerSize', 8,'MarkerEdgeColor','k', 'MarkerFaceColor', cmap(color_ind*60,:)); hold on

            color_ind = color_ind - 1; % coincidence with the vector plot

        end
        xlim([-5 -4]); ylim([-0.4, 0])

        savename = strcat('y_0=', title_y0, '_theta_0=', title_deg, '_orientation_vs_x.eps');

        set_plot(gcf, gca) 
%         f=gcf;
%         exportgraphics(f,['F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\' ...
%             'Figures\about orientations vs x\until contact\Given y0 and theta0\', savename{1}]);

%         close
    end
end


