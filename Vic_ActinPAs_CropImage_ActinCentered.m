%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Center the filament in cropped image for AI training.

%  Good_case_frm: stored the index of the 'Good_case' variable. Get from
%                 the 'draw reconstruction and selection' process.
%                 * call in loop j: xy(1).frame(Good_case_frm(j))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

crop_size = 512;
% save_pathname = ['F:\Experimental Data (EXTRACTED)\Actin Filaments in ' ...
%     'Porous Media\Cropped images for AI tracking\'];
save_pathname_uint8 = ['F:\Experimental Data (EXTRACTED)\Actin Filaments in ' ...
    'Porous Media\Cropped images for AI tracking (uint8)\'];

xlsfile = readcell('ForActinPostprocessing.xlsx','Sheet','Sheet1','NumHeaderLines',1);
% This is the file that contains all the information about the later processing (in sheet 1).

ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.

for no_Group = [7 8 13:28]

    the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
    thefiles = dir(fullfile(storePath{no_Group},'*.mat'));

    theimgs = dir(['Z:\Experimental Data (EXTRACTED)\Actin Filaments in Porous' ...
        ' Media\',num2str(the_exp_date),'-Actin\*.tif']);

    for file_ind = 1:length(thefiles)

        filename = thefiles(file_ind).name;

        if contains(filename, 'PAsInfoAdded_')

            load(fullfile(thefiles(1).folder, thefiles(file_ind).name));

            filename = thefiles(file_ind).name
            save_filename = [num2str(the_exp_date), filename(38:end-17), '_frm'];

            % find the cooresponding *.tif image
            image_names = struct2cell(theimgs); image_names = image_names(1, :);
            image_ind = find(cellfun(@(x) contains(x, filename(25:end-17)), image_names));
            if length(image_ind) ~= 1
                image_ind = find(cellfun(@(x) strcmp(x(1:end-4), filename(25:end-17)), image_names));
            end

            for frm_ind = 1:size(Good_case_frm,2)

                xy_ind = Good_case_frm(frm_ind); % index of the 'good' cases
                II = imread(fullfile(theimgs(1).folder, theimgs(image_ind).name), ...
                    xy(1).frame(xy_ind)); % load the image

                CoM_xy = xy.centroid{1,xy_ind}; CoM_xy(2) = 2048-CoM_xy(2);
                if no_Group == 25
                    CoM_xy(1) = 2048-CoM_xy(1);
                end
                if round(CoM_xy(1)-crop_size/2)>0 && round(CoM_xy(2)-crop_size/2)>0 && round(CoM_xy(1)+crop_size/2)<2049 && round(CoM_xy(2)+crop_size/2)<2049
                    Cropped_II = II(round(CoM_xy(2)-crop_size/2):round(CoM_xy(2)+crop_size/2)-1, ...
                        round(CoM_xy(1)-crop_size/2):round(CoM_xy(1)+crop_size/2)-1);
                    if ~isa(Cropped_II,'uint8')
                        Cropped_II_uint8 = uint8(double(Cropped_II)/double(max(max(Cropped_II))) * 255); 
                        % Convert to 8-bits by linearly scaling from min-max to 0-255 (similar to ImageJ)
                    else
                        Cropped_II_uint8 = Cropped_II;
                    end

%                     imwrite(Cropped_II, [save_pathname_uint16, save_filename, num2str(xy_ind), '.tif']);
                    imwrite(Cropped_II_uint8, [save_pathname_uint8, save_filename, num2str(xy_ind), '.tif']);
                end

            end

            clearvars CoM_xy Cropped_II Cropped_II_uint8
        end
    end
end
