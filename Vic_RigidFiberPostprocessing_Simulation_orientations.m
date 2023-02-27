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
            ylabel('$Orientation\ (^o)$','FontSize', 18,'Interpreter', 'latex');

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


%%  
clear; close all; clc;

Files = dir(['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared' ...
    '\Data_Give_to_Zhibo_20230223_videos\Orientations\*.mat']);

n = 1;
for ii = 1: length(Files)

    load(fullfile(Files(ii).folder, Files(ii).name));
    current_deg =  extractBetween(Files(ii).name,'theta0_', '_L');
    current_deg = strrep(current_deg,'o','.'); current_deg = strrep(current_deg,'m','-');
    current_L =  extractBetween(Files(ii).name,'L_', '_y0');
    current_L = strrep(current_L,'o','.');
    current_y0 =  extractBetween(Files(ii).name,'y0_', '_ori');
    current_y0 = strrep(current_y0,'o','.');
    if isnan(str2double(current_y0))
        current_y0 = current_y0(1:end-10);
    end

    if str2double(current_deg) ~= 2.5 || str2double(current_L) ~= 0.5
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

figure('color', 'w'); set(gcf, 'Position', [100 100 1000 800]);
cmap = cmocean('thermal');  legend_txt = {};
[~, sortID] = sort(cell2mat(To_Plot(:,3)));  % sort the plotting order

color_ind = 1;
for jj = [3 5 6 7 8 9 13 16]%1:size(To_Plot, 1)

    plot(To_Plot{sortID(jj), 2}, To_Plot{sortID(jj), 1}, 'o','MarkerSize', 5,'MarkerEdgeColor','k', ...
        'MarkerFaceColor', cmap(color_ind*30,:)); hold on

    xlabel('$x/h_{obs}$','FontSize', 18,'Interpreter', 'latex');
    ylabel('$Orientation\ (^o)$','FontSize', 18,'Interpreter', 'latex');

    title(strcat('$\theta_0=', title_deg, '\ and\ L=', title_L, '$'), ...
        'FontSize', 24, 'Interpreter', 'latex');
    legend_txt = [legend_txt, strcat('$y_0=', num2str(To_Plot{sortID(jj), 3}),'$')];

    color_ind = color_ind + 1;

end
xlim([-10 10])
ax = gca; ax.FontSize = 18;
legend(legend_txt, 'Location', 'northwest', 'FontSize', 18,  'Interpreter', 'latex');
