%%%% Calculate contour length, mu_bar and mean(U_x)/U_0 %%%%
% data from Vic_ActinAddInformationSave.m
% data name format: PAsInfoAdded_trajectory_..._.mat
% results save in the Excels named: Dynamics xxxxxxxx-Actin.xlsx

%%%% The contour L is the average of 10% longest snapshots except for the
%%%% extreme values (out of 2-sigma among those averaged ones). 

clear; close all; clc;
f = figuresetting('centimeters',21,15,'Helvetica',24,'off',2,'off','off');
mag = 0.1; % um/pixel

xlsfile = readcell('ForActinPostprocessing.xlsx','Sheet','Sheet1','NumHeaderLines',1);
% This is the file that contains all the information about the later processing (in sheet 1).

NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.
Init_U = xlsfile(:, 8);  % initial velocity, unit: m/s.
D_pillar = xlsfile(:, 13);  % pillar diameter, unit: um

excelpathname = ['F:\Processing & Results\Actin Filaments in Porous Media' ...
    '\Dynamics manually classification\'];
% The path contains the Excels be be written in.

for no_Group = [7 8 13:28]

    the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
    thefiles = dir(fullfile(storePath{no_Group},'*.mat'));

    excelname = ['Dynamics ', num2str(the_exp_date), '-Actin.xlsx'];
    xlsfile_to_write = readcell([excelpathname, excelname],'Sheet','Sheet1','NumHeaderLines',1); 

    all_names = xlsfile_to_write(:, 1);

    for file_ind = 1:length(thefiles)

        filename = thefiles(file_ind).name;

        if contains(filename, 'PAsInfoAdded_') 

            load(fullfile(thefiles(1).folder, thefiles(file_ind).name));

            pathname = thefiles(1).folder;
            filename = thefiles(file_ind).name;

            excel_pos = find(cellfun(@(x) contains(filename, x), all_names));

            ContourL_all = xy.arclen_spl(Good_case_frm);
            ContourL = VicFc_Get_ContourLength(ContourL_all) * mag; % unit: um

            mu_bar = VicFc_Get_elastoviscousNum(ContourL*1e-6, Init_U{no_Group}, D_pillar{no_Group}*1e-6);

            if contains(filename, 'interval2')
                Ux_to_U0 = (xy.centroid{1, Good_case_frm(end)}(1) - ...
                    xy.centroid{1, Good_case_frm(1)}(1)) * mag * 1e-6 / ...
                ((xy(1).frame(Good_case_frm(end)) - xy(1).frame(Good_case_frm(1))) * 0.04) / Init_U{no_Group};
            else
                Ux_to_U0 = (xy.centroid{1, Good_case_frm(end)}(1) - ...
                    xy.centroid{1, Good_case_frm(1)}(1)) * mag * 1e-6 / ...
                ((xy(1).frame(Good_case_frm(end)) - xy(1).frame(Good_case_frm(1))) * 0.02) / Init_U{no_Group};
            end
            
            Loc = ['B', num2str(excel_pos+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
            writematrix(ContourL,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value inti the excel.

            Loc = ['N', num2str(excel_pos+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
            writematrix(mu_bar,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value inti the excel.

            Loc = ['S', num2str(excel_pos+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
            writematrix(Ux_to_U0,[excelpathname, excelname],'Sheet','Sheet1','Range', Loc);  % Write the value inti the excel.

%             % to plot the distribution of the lengths
%             f.figure('');
%             h = histogram(ContourL,round(numel(ContourL)/10));
%             f.interp_font('latex');
%             f.axes('linear',[min(ContourL) max(ContourL)],'linear',[0 4*round(numel(ContourL)/10)],...
%                 ' $L$ $(\mu m)$ ',' $Case\ number$ ',24);
% 
%             % to plot the max, min, and mean of the lengths 
%             f.figure('');
%             boxplot(ContourL);
%             f.interp_font('latex');
%             f.axes('linear',[0 2],'linear',[min(ContourL)-50 max(ContourL)+50],...
%                 ' ',' $L$ $(\mu m)$ ',24);
% 
%             close all
            
            clearvars xy Good_case_frm
        end 
    end
end
