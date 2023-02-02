%%%% Add information to the trajectory_M63_..._batch1.mat of filament.
% 1. lzero
% 2. Good_case_frm
% 3. pillar array locations
%%%% because the information is different from one case to another.

clear; clc; close all;
xlsfile = readcell('ForActinPostprocessing.xlsx','Sheet','Sheet1','NumHeaderLines',1); % This is the file that contains all the information about the later processing (in sheet 1).

NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
storePath = xlsfile(:, 2);  % Path of the data to be processed.
PAsPath = xlsfile(:, 3);  % Path of the pillar array information.

for no_Group = 1: NumGroup

    filelist = dir(fullfile(storePath{no_Group},'*.mat'));
    % list of the .mat files which contain the reconstruction information
    % (came from 'Filaments detection' code) in one group.

    if no_Group <= 12 || no_Group == 16
        % experiments before 2022.01.04 and experiments on 2022.02.17 Group1

        for no_Case = 1:length(filelist)

            load(PAsPath{no_Group}); % Load the pillar array information.
            load([storePath{no_Group}, filesep , filelist(no_Case).name])
            filename = filelist(no_Case).name;
            save_filename = ['PAsInfoAdded_', filename];

            Good_case_frm = Good_case;
            lzero = zeros(size(Good_case_frm,2), 1);
            for frm_ind = 1:size(Good_case_frm,2)
                xy_ind = Good_case_frm(frm_ind);
                lzero(frm_ind) = max(lobject,ceil(5*lnoise));
            end

            % save ...
%             save([storePath{no_Group}, filesep , save_filename], 'ds','FilNum', ...
%                 'final_frame','frame_step','framelist','Good_case','improc', ...
%                 'InfoImage','initial_frame','lnoise','lobject', 'MinBranchLength', ...
%                 'missed_frames','N_fil','npnts','prcs_img','ROI','sensitivity', ...
%                 'structsensitivity','thickness','threshold','xskip','yskip', ...
%                 'xy','xwin','ywin','centers','circleMask','metric','radii','Good_case_frm','lzero')

            % clear the saved variables
            clearvars -except PAsPath storePath NumGroup no_Group no_Case filelist
        end

    elseif no_Group == 13 || no_Group == 14 || no_Group == 15 || no_Group == 17 || no_Group == 18
        % experiments on 2022.02.16 Group1, Group2 & Group3 and experiments on 2022.02.17 Group2 & Group3

        for no_Case = 1:length(filelist)

            load(PAsPath{no_Group}); % Load the pillar array information.
            load([storePath{no_Group}, filesep , filelist(no_Case).name])
            filename = filelist(no_Case).name;
            save_filename = ['PAsInfoAdded_', filename];

            Good_case_frm = find(ismember(xy(1).frame, Good_case));
            lzero = zeros(size(Good_case_frm,2), 1);

            % save ...
%             save([storePath{no_Group}, filesep , save_filename], 'framelist', ...
%                 'Good_case','InfoImage','prcs_img','prmt','ROI','xy','centers', ...
%                 "circleMask",'metric','radii','Good_case_frm','lzero')

            % clear the saved variables
            clearvars -except PAsPath storePath NumGroup no_Group no_Case filelist
        end
    end
end








%% Draw and check if the information added correctly.
clear; close all; clc;
xlsfile = readcell('ForActinPostprocessing.xlsx','Sheet','Sheet1','NumHeaderLines',1); 
% This is the file that contains all the information about the later processing (in sheet 1).

NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.
PAsPath = xlsfile(:, 3);  % Path of the pillar array information.

for no_Group = 1: NumGroup

    the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
    filelist = dir(fullfile(storePath{no_Group},'*.mat'));
    % list of the .mat files which contain the reconstruction information
    % (came from 'Filaments detection' code) in one group.

    theimgs = dir(['F:\Experimental Data (EXTRACTED)\Actin Filaments in Porous' ...
        ' Media\',num2str(the_exp_date),'-Actin\AfterAveBGR\*.tif']);
    % The *.tif files where stores your results.
    
    for no_Case = 1:length(filelist)

        filename = filelist(no_Case).name;

        if contains(filename, 'PAsInfoAdded_')
            load([storePath{no_Group}, filesep , filename])

            Foo = ceil(length(Good_case_frm)/2);
            xy_ind = Good_case_frm(Foo);% index of the 'good' cases

            image_names = struct2cell(theimgs); image_names = image_names(1, :);
            image_ind = find(cellfun(@(x) contains(x, filename(25:end-11)), image_names));
            II = imread(fullfile(theimgs(1).folder, theimgs(image_ind).name), xy(1).frame(xy_ind)); % load the image

            fiber_center = xy.centroid{xy_ind}; 
            fiber_center_col = round(fiber_center(1) + lzero(Foo));
            fiber_center_row = round(2048 - fiber_center(2) - lzero(Foo));
            row_up = max(fiber_center_row-200, 1); row_down = min(fiber_center_row+200, 2048);
            col_left = max(fiber_center_col-200, 1); col_right = min(fiber_center_col+200, 2048);
            imshow(II(row_up:row_down, col_left:col_right), []); hold on
            plot(xy.spl{xy_ind}(:,1)+lzero(Foo)-col_left+1, ...
                2048-xy.spl{xy_ind}(:,2)-lzero(Foo)-row_up+1, 'LineWidth', 2); 

            f=gcf;
            exportgraphics(f,['D:\Dropbox\GitHub\tmp', filesep, num2str(the_exp_date), filename(25:end-11), '.png'],'Resolution',100)
            close
        end       
    end
end



