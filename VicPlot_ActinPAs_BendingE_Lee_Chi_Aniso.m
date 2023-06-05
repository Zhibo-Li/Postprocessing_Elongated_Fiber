%%%% Plot the chronophotograoh, bending energy and Lee of actin filaments.
% data from Vic_ActinPAs_BendingE_Lee.m
% data name format: PlusInfo_..._.mat

clear; close all; clc;

set(0,'DefaultAxesFontSize',14);
set(0,'defaulttextfontsize',14);
% set(0,'defaultAxesFontName', 'times new roman');
% set(0,'defaultTextFontName', 'times new roman');
set(0,'defaultAxesFontName', 'Arial');
set(0,'defaultTextFontName', 'Arial');
set(0,'defaulttextInterpreter','latex') 

mag = 0.1; % um/pixel

xlsfile = readcell('ForActinPostprocessing.xlsx','Sheet','Sheet1','NumHeaderLines',1);
% This is the file that contains all the information about the later processing (in sheet 1).

NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.

for no_Group = [7 8 13:28]

    the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
    theTruestorePath = storePath{no_Group};
    theTruestorePath = strrep(theTruestorePath,'results','results_plus');
    thefiles = dir(fullfile(theTruestorePath,'*.mat'));

    for file_ind = 1:length(thefiles)

        filename = thefiles(file_ind).name;

        if contains(filename, 'PlusInfo_') 

            load(fullfile(thefiles(1).folder, thefiles(file_ind).name));

            pathname = thefiles(1).folder;
            filename = thefiles(file_ind).name

            figure('color', 'w'); set(gcf, 'Position', [100 300 1200 400]);

            yyaxis left
            plot((CoM_x+lzero)*mag, Energy, '.b', 'LineStyle','none', 'MarkerSize', 20);
            ylabel('Bending energy', 'FontName', 'Times new roman');

            yyaxis right         
            plot((CoM_x+lzero)*mag, L_ee_norm,'dr', 'LineStyle','none', 'MarkerSize', 8);
            hold on
            plot((CoM_x+lzero)*mag, Chi/pi,'sm', 'LineStyle','none', 'MarkerSize', 8);
            hold on
            plot((CoM_x+lzero)*mag, aniso, '^k', 'LineStyle','none', 'MarkerSize', 8); 
            ylim([-0.5 1.5]);
            ylabel({'End-to-end distance, ','orientation, and sphericity'}, 'FontName', 'Times new roman');

            legend('$E_{Bending}\ (J)$', '$L_{ee}/L_0$', '$\chi / \pi$', ...
                '$\omega$', 'FontSize', 18,'Interpreter', 'latex','NumColumns', 4, ...
                'location', 'northwest');
%           Sphericity:  \omega = 1 - 4{\lambda}_1{\lambda}_2/({\lambda}_1+{\lambda}_2)^2

            xlim([0 2050]*mag);
            xlabel('$x\ (\mu{m})$')
            grid on; hold on

            savepath = ['F:\Processing & Results\Actin Filaments in Porous Media\Figures\Energy_Lee_Chi_Aniso',...
                pathname(56:70), pathname(end-7:end)]; 
            if ~exist(savepath, 'dir')
                mkdir(savepath)
            end

            set(gca, 'Box', 'On', 'FontSize', 18, 'TickLabelInterpreter','latex')

            f=gcf;
            exportgraphics(f,[savepath, filesep, filename(1: end-4), '_E_Lee_Chi_Aniso.png'],'Resolution',100)

            close all
            clearvars CoM_x Energy xy centers radii L_ee_norm Chi aniso
        end 
    end
end







