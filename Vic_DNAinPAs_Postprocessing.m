%% Add the stage positions from tracking.

clear; clc; close all;

objective_magnification = 0.16; % the magnification of the objective (um/pixel)

% select the folder where the data is stored.
dataPath = uigetdir('Z:\Processing & Results\DNA in Pillar Arrays', ...
    'Select the folder where the data is stored.');
stageInfoPath = uigetdir('Z:\Experimental Data (RAW)\DNA in Pillar Arrays', ...
    'Select the folder where the stage information is stored.');
figurePath = uigetdir('Z:\Processing & Results\DNA in Pillar Arrays', ...
    'Select the folder where the figures to be saved.');
filelist = dir(fullfile(dataPath,'*.mat'));
% list of the .mat files which contain the reconstruction information

% load the *.mat file in a loop
for Test_No = 1:length(filelist)
    
    filename = filelist(Test_No).name;
    if startsWith(filename, 'trajectory_')
        
        %         % The absolute positions after considering the stage movement.
        %         save_filename = ['AbsPosition_', filename];
        
        % load the file
        load([dataPath, filesep , filename]);
        
        % extract the char between 'trajectory_' and '_batch1.mat' from filename
        case_name = extractBetween(filename, 'trajectory_', '_batch1.mat');
        
        % The corresponding RAW tifs
        rawdata_tifs = dir(fullfile(stageInfoPath, filesep, case_name{1},'*.tif'));
        img_no_str = cat(1, rawdata_tifs(:).name);
        img_no = str2num(img_no_str(:, 6:12)); % extract the image number
        
        % load the stage information (*.txt file)
        StagePosition_txt = readmatrix([stageInfoPath, filesep, case_name{1}, filesep, 'StagePosition.txt'],'Range','A:D');
        StagePosition_txt(~ismember(StagePosition_txt(:,1), img_no), :) = []; % remove the invalid rows
        if img_no(1) ~= StagePosition_txt(1,1)
            StagePosition_txt = [StagePosition_txt(1, :); StagePosition_txt];
            StagePosition_txt(1) = StagePosition_txt(1) - 1; % add the missed FIRST frame
        end
        if ~isempty(img_no(~ismember(img_no, StagePosition_txt(:,1))))
            missed_frame = find(img_no == img_no(~ismember(img_no, StagePosition_txt(:,1)))); % find the missed frame
        else
            missed_frame = [];
        end
        if length(img_no) < size(StagePosition_txt, 1) + length(missed_frame)
            StagePosition_txt(diff(StagePosition_txt(:,1))==0,:) = []; % remove the duplicate
        end
        
        Good_case(ismember(Good_case, missed_frame))= []; % no information about the missed frame
        
        false_ind = ~ismember(xy(1).frame, Good_case);
        % Here, the 'Good_case' means the absolute frame number, not the index of
        % xy.frame (different to the code before which calculation and selection
        % are two processes.
        
        xy.crd(false_ind) = [];
        xy.centroid(false_ind) = [];
        xy.arclen(false_ind) = [];
        xy.seglen(false_ind) = [];
        xy.nGoodframe = length(Good_case) - xy(1).nframe + xy(1).nframe;
        xy.frame(false_ind) = [];
        xy.spl(false_ind) = [];
        xy.knots(false_ind) = [];
        xy.arclen_spl(false_ind) = [];
        xy.seglen_spl(false_ind) = [];
        
        img_no_Good_case = img_no(Good_case, :);
        StagePosition_Good_case = StagePosition_txt;
        StagePosition_Good_case(~ismember(StagePosition_Good_case(:,1), ...
            img_no_Good_case), :) = [];
        
        StageXY_Good_case = StagePosition_Good_case(:, 2:3);
        StageXY_relative_Good_case = StageXY_Good_case - StageXY_Good_case(1, :);
        
        for Good_case_No = 1:length(Good_case)
            
            xy.crd{Good_case_No} = xy.crd{Good_case_No} + StageXY_relative_Good_case(Good_case_No, :)/objective_magnification;
            xy.centroid{Good_case_No} = xy.centroid{Good_case_No} + StageXY_relative_Good_case(Good_case_No, :)/objective_magnification;
            xy.spl{Good_case_No} = xy.spl{Good_case_No} + StageXY_relative_Good_case(Good_case_No, :)/objective_magnification;
            xy.knots{Good_case_No} = xy.knots{Good_case_No} + StageXY_relative_Good_case(Good_case_No, :)/objective_magnification;
            
        end
        
        % save ...
        %                 save([storePath{no_Group}, filesep , save_filename], 'framelist', ...
        %                     'Good_case','InfoImage','prcs_img','prmt','ROI','xy','centers', ...
        %                     "circleMask",'metric','radii','Good_case_frm','lzero')
        
        %         % clear the saved variables
        %         clearvars -except PAsPath storePath NumGroup no_Group Test_No filelist
    end
    

    figure('color', 'w'); set(gcf, 'Position', [100 300 1600 100]);
    for k = 1:length(Good_case)
        plot(xy.spl{k}(:,1)*objective_magnification, xy.spl{k}(:,2)*objective_magnification); hold on
    end
    xlim([0 max(xy.spl{k}(:,1)*objective_magnification)]); ylim([0 150])
    xlabel('$x\ (\mu{m})$', 'Interpreter', 'latex'); ylabel('$y\ (\mu{m})$', 'Interpreter', 'latex');

    exportgraphics(gcf,[figurePath, filesep, case_name{1}, '_Chronophotography.png'],'Resolution',100)
    

    figure('color', 'w'); set(gcf, 'Position', [100 300 1500 200]);
    V_front_group = zeros(1, length(Good_case)-1);
    centroid_x_group = zeros(1, length(Good_case)-1);
    for k = 1:length(Good_case)
        
        centroid_x = xy.centroid{k}(:,1)*objective_magnification;
        
        % calculate the end to end distance of the elongated fiber
        L_ee = norm(xy.spl{k}(1,:)-xy.spl{k}(end,:))*objective_magnification;
        
        % calculate the distances between the centroid and the ends of the fiber
        D_end_1 = norm(xy.spl{k}(1,:)-xy.centroid{k})*objective_magnification;
        D_end_2 = norm(xy.spl{k}(end,:)-xy.centroid{k})*objective_magnification;
        if xy.spl{k}(1,1) > xy.spl{k}(end,1)
            D_front = D_end_1; D_rear = -D_end_2;
            if k > 1
                V_front_group(k-1) = norm(xy.spl{k}(1,:)-xy.spl{k-1}(1,:))*objective_magnification/(img_no_Good_case(k)-img_no_Good_case(k-1));
                centroid_x_group(k-1) = centroid_x;
            end
        else
            D_front = D_end_2; D_rear = -D_end_1;
            if k > 1
                V_front_group(k-1) = norm(xy.spl{k}(end,:)-xy.spl{k-1}(end,:))*objective_magnification/(img_no_Good_case(k)-img_no_Good_case(k-1));
                centroid_x_group(k-1) = centroid_x;
            end
        end
        
        yyaxis left
        plot(centroid_x, D_front, 'Color','m', 'LineStyle','none', ...
            'Marker','*', 'MarkerSize',3); hold on
        plot(centroid_x, D_rear, 'Color','b', 'LineStyle','none', ...
            'Marker','^', 'MarkerSize',3); hold on
        ylabel('$d_{\rm end-to-center}\,(\mu{m})$', 'Interpreter', 'latex');
        ylim([-15 15])
        
        yyaxis right
        plot(centroid_x, L_ee, 'Color','k', 'LineStyle','none', ...
            'Marker','o', 'MarkerSize',3); hold on
        ylabel('$L_{\rm ee}\,(\mu{m})$', 'Interpreter', 'latex');
        ylim([0 30])
    end
    
    % plot the position of the pillars
    pillar_start = 435*objective_magnification;
    pillar_end = max(centroid_x);
    pillar_interval = 30;
    for pillar_x = pillar_start:pillar_interval:pillar_end
        hold on;
        line([pillar_x pillar_x], [0 150], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 0.5);
    end
    %% 
    
    xlim('auto');
    legend({'Front end','Rear end'}, 'FontName', 'Times New Roman');
    xlabel('$x\,(\mu{m})$', 'Interpreter', 'latex');
    
    exportgraphics(gcf,[figurePath, filesep, case_name{1}, '_Lee.png'],'Resolution',100)

    
    figure('color', 'w'); set(gcf, 'Position', [100 300 1500 200]);
    plot(centroid_x_group, V_front_group, 'Color','k', 'LineStyle','none', ...
        'Marker','.', 'MarkerSize',10); hold on
    xlabel('$x\,(\mu{m})$', 'Interpreter', 'latex');
    ylabel('$V_{\rm front}\,(a.u.)$', 'Interpreter', 'latex');
    ylim([0 15])
    for pillar_x = pillar_start:pillar_interval:pillar_end
        hold on;
        line([pillar_x pillar_x], [0 150], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 0.5);
    end
    
    exportgraphics(gcf,[figurePath, filesep, case_name{1}, '_V_front.png'],'Resolution',100)
end
