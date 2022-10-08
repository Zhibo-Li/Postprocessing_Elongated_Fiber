%% Add 'real time' to the results
% * For the rigid fiber *
% The time data is extracted from the *.cxd file and saved to the *.xlsx

clear; close all; clc;

% input the table of time.
[filename, pathname] = uigetfile(['F:\Experimental Data (EXTRACTED)\FSI - Ri' ...
    'gid Fiber &  Individual Obstacle\*.xlsx'], 'Choose the EXCEL file');
thetime = readmatrix(fullfile(pathname, filename), "NumHeaderLines", 1);
timecase_No = size(thetime, 2);
thefilenames = readcell(fullfile(pathname, filename), "Range", "A1:ZZ1");
thefilenames = thefilenames(1, 1:timecase_No); 

% input the *.mat file .
matpath = uigetdir('F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\', ...
    'Choose the folder containing the *.mat result files');
matFiles = dir(fullfile(matpath, '*.mat'));

% folder to save the results
savepathname = uigetdir(['F:\Processing & Results\FSI - Rigid Fiber &  Individua' ...
    'l Obstacle\'], 'Choose a folder to save the results');

% loop the excel columns to add the git a.
for ii = 1:timecase_No 

    tmp_currenttime = thetime(:, ii);
    tmp_currenttime(isnan(tmp_currenttime)) = [];
    timestamps = tmp_currenttime - min(tmp_currenttime);  % the timestamps of a given case.
     
    matcase_No = find(contains({matFiles.name}, thefilenames{ii}));
    %%% !!! This condition is not enough, e.g., AAA_3 and AAA_32 would be selected together.
    if isempty(matcase_No)
        continue
    end

    for jj = 1:length(matcase_No)

        currenttime = timestamps;

        load([matFiles(matcase_No(jj)).folder, filesep, matFiles(matcase_No(jj)).name]);
        currentmatfile = matFiles(matcase_No(jj)).name;  % find and load the counterpart *.mat file.

        if contains(currentmatfile, '_no') && contains(currentmatfile, '-no')
            if isnan(str2double(extractBetween(currentmatfile,'_no','-no')))
                imageinterval = str2double(extractBetween(currentmatfile,'-interval','-no'));
                imagestart = str2double(extractBetween(currentmatfile,'_no','-interval'));
                currenttime(1:imagestart-1) = []; % remove the non-calculated counterpart timestamps.
                currenttime(2:imageinterval:end) = [];
            else
                imagestart = str2double(extractBetween(currentmatfile,'_no','-no')); % the location of the start frame in the original cxd file.
                currenttime(1:imagestart-1) = []; % remove the non-calculated counterpart timestamps.
            end
        end

        calculated_frame_time = currenttime(xy.frame);
        Good_case_frm_time = calculated_frame_time(Good_case_frm); % add the timestamps

        if ~exist(savepathname,'dir')
            mkdir(savepathname)
        end

        save(fullfile(savepathname, currentmatfile),'Good_case_frm_time','Good_case_frm','xy','ROI','prmt');
    end

    clearvars matcase_No

end