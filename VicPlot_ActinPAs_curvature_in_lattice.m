%% Plot fiber curvatures in the lattice.
% data from Vic_ActinPAs_Normalize_fiber_xy_in_lattice.m
% data name format: PlusInfo_trajectory_..._batch1.mat

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
Array_angles = cell2mat(xlsfile(:, 14));  % The flow angles.

plot_resolution = 50;
Xedges = [0:1/plot_resolution:1]; Yedges = [0:1/plot_resolution:1]; % for the histogram

xlsFolder = dir(['F:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Dynamics manually classification\*.xlsx']);
% Folder contains 'Dynamics manually classification' Excels (to remove bad reconstructions)

for current_angle = [0 10 15 20 30 35 45] % different angles

    fiber_xy_in_lattice = [];
    curvature = [];

    if current_angle == 0
        Groups = [7;8];
    elseif current_angle == 20
        Groups = [16;17;18];
    else
        Groups = find(Array_angles == current_angle);
    end

    for ii = 1: length(Groups) % different exp. groups (might be different days)

        no_Group = Groups(ii);

        the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
        theTruestorePath = storePath{no_Group};
        theTruestorePath = strrep(theTruestorePath,'results','results_plus');
        thefiles = dir(fullfile(theTruestorePath,'*.mat'));

        Dyn_xlsfile = readcell([xlsFolder(1).folder, filesep, 'Dynamics ', num2str(the_exp_date),'-Actin.xlsx'], ...
        'Sheet','Sheet1','NumHeaderLines',1);
        % load 'Dynamics manually classification' Excel
        case_names = Dyn_xlsfile(:, 1);
        if_has_note = Dyn_xlsfile(:, 11);
        if_out_of_plane = Dyn_xlsfile(:, 12);
        if_bad_reconstruction = Dyn_xlsfile(:, 13);
        BadCaseInd = cell2mat(if_has_note) | cell2mat(if_out_of_plane) | cell2mat(if_bad_reconstruction); 
        % The index of bad reconstructed cases.
        bad_case_names = case_names(BadCaseInd);
        % The names of bad reconstructed cases.

        for file_ind = 1:length(thefiles) % load the cases

            filename = thefiles(file_ind).name;

            if ~any(strcmp(bad_case_names, filename(10: end-4))) % Remove bad reconstructions
                if contains(filename, 'PlusInfo_')

                    load(fullfile(thefiles(1).folder, thefiles(file_ind).name));

                    for frm_ind = 1:size(Good_case_frm,2) % different snapshots

                        xy_ind = Good_case_frm(frm_ind); % index of the 'good' cases
                        fiber_xy_in_lattice = [fiber_xy_in_lattice; xy.fiber_xy_in_lattice{1, xy_ind}];
                        curvature = [curvature; xy.curvature{1, xy_ind}/mag]; % /mag: convert curvature unit to 1/um

                    end
                end
            end
        end
    end

    fiber_xy_in_lattice(fiber_xy_in_lattice > 1) = 1;
    fiber_xy_in_lattice(fiber_xy_in_lattice < 0) = 0;  % to make sure no error in accumarray
    [~,~,~,ind_x,ind_y] = histcounts2(fiber_xy_in_lattice(:,1), ...
        fiber_xy_in_lattice(:,2), Xedges, Yedges);
    binned = accumarray([ind_y,ind_x],curvature,[numel(Xedges)-1 numel(Xedges)-1],@VicMean,0);

    [x_plot,y_plot] = meshgrid(0.5/plot_resolution:1/plot_resolution:1-0.5/plot_resolution, ...
        0.5/plot_resolution:1/plot_resolution:1-0.5/plot_resolution);
    figure('color', 'w'); set(gcf, 'Position', [100 100 600 600]);
    
    pcolor(x_plot,y_plot, binned);
    axis equal; axis off
    cmocean('deep'); caxis([0 1]); %colorbar;
    
    viscircles([0 0; 0 1; 1 0; 1 1], mean(radii)*ones(4, 1)/mean([Ctr2Ctr_x, Ctr2Ctr_y]),'Color','r');
    xlim([0 1]); ylim([0 1]);

    f = gcf;
    savefig(f,['F:\Processing & Results\' ...
        'Actin Filaments in Porous Media\Figures\Curvature Distribution in Lattice\' ...
        'actin_curvature_lattice_angle_',num2str(current_angle),'_noBadCases.fig'])

    set(f,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
        'Actin Filaments in Porous Media\Figures\Curvature Distribution in Lattice\' ...
        'actin_curvature_lattice_angle_',num2str(current_angle),'_noBadCases.eps']);

    clearvars fiber_xy_in_lattice curvature

end

% another 20 degree:
no_Group = 28;
fiber_xy_in_lattice = [];
curvature = [];

the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
theTruestorePath = storePath{no_Group};
theTruestorePath = strrep(theTruestorePath,'results','results_plus');
thefiles = dir(fullfile(theTruestorePath,'*.mat'));

Dyn_xlsfile = readcell([xlsFolder(1).folder, filesep, 'Dynamics ', num2str(the_exp_date),'-Actin.xlsx'], ...
    'Sheet','Sheet1','NumHeaderLines',1);
% load 'Dynamics manually classification' Excel
case_names = Dyn_xlsfile(:, 1);
if_has_note = Dyn_xlsfile(:, 11);
if_out_of_plane = Dyn_xlsfile(:, 12);
if_bad_reconstruction = Dyn_xlsfile(:, 13);
BadCaseInd = cell2mat(if_has_note) | cell2mat(if_out_of_plane) | cell2mat(if_bad_reconstruction);
% The index of bad reconstructed cases.
bad_case_names = case_names(BadCaseInd);
% The names of bad reconstructed cases.

for file_ind = 1:length(thefiles) % load the cases

    filename = thefiles(file_ind).name;

    if ~any(strcmp(bad_case_names, filename(10: end-4))) % Remove bad reconstructions
        if contains(filename, 'PlusInfo_')
            load(fullfile(thefiles(1).folder, thefiles(file_ind).name));

            for frm_ind = 1:size(Good_case_frm,2) % different snapshots

                xy_ind = Good_case_frm(frm_ind); % index of the 'good' cases
                fiber_xy_in_lattice = [fiber_xy_in_lattice; xy.fiber_xy_in_lattice{1, xy_ind}];
                curvature = [curvature; xy.curvature{1, xy_ind}/mag]; % /mag: convert curvature unit to 1/um

            end
        end
    end
end

fiber_xy_in_lattice(fiber_xy_in_lattice > 1) = 1;
fiber_xy_in_lattice(fiber_xy_in_lattice < 0) = 0;  % to make sure no error in accumarray
[~,~,~,ind_x,ind_y] = histcounts2(fiber_xy_in_lattice(:,1), ...
    fiber_xy_in_lattice(:,2), Xedges, Yedges);
binned = accumarray([ind_y,ind_x],curvature,[numel(Xedges)-1 numel(Xedges)-1],@VicMean,0);

[x_plot,y_plot] = meshgrid(0.5/plot_resolution:1/plot_resolution:1-0.5/plot_resolution, ...
    0.5/plot_resolution:1/plot_resolution:1-0.5/plot_resolution);
figure('color', 'w'); set(gcf, 'Position', [100 100 600 600]);

pcolor(x_plot,y_plot, binned);
axis equal; axis off
cmocean('deep'); caxis([0 1]); %colorbar;

viscircles([0 0; 0 1; 1 0; 1 1], mean(radii)*ones(4, 1)/mean([Ctr2Ctr_x, Ctr2Ctr_y]),'Color','r');
xlim([0 1]); ylim([0 1]);

f = gcf;
savefig(f,['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Curvature Distribution in Lattice\' ...
    'actin_curvature_lattice_angle_20_smallerpillar_noBadCases.fig'])

set(f,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Curvature Distribution in Lattice\' ...
    'actin_curvature_lattice_angle_20_smallerpillar_noBadCases.eps']);



%% Plot only the colorbar for the above figures.
figure('color', 'w'); set(gcf, 'Position', [100 100 200 380]);

cmocean('deep');
caxis([0 1]);
c = colorbar;
c.Label.String = '$\rm{Curvature}\, (1/\mu m)$';
c.Label.Interpreter = 'LaTeX';
c.TickLabelInterpreter = 'LaTeX';
c.FontSize = 20;
c.Location = 'west';
axis off

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
        'Actin Filaments in Porous Media\Figures\Curvature Distribution in Lattice\' ...
        'actin_curvature_lattice_ColorBar.eps']);



%% Plot fiber curvatures in the lattice (L < 30um).
% data from Vic_ActinPAs_Normalize_fiber_xy_in_lattice.m
% data name format: PlusInfo_trajectory_..._batch1.mat

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
Array_angles = cell2mat(xlsfile(:, 14));  % The flow angles.

plot_resolution = 50;
Xedges = [0:1/plot_resolution:1]; Yedges = [0:1/plot_resolution:1]; % for the histogram

for current_angle = [0 10 15 20 30 35 45] % different angles

    fiber_xy_in_lattice = [];
    curvature = [];

    if current_angle == 0
        Groups = [7;8];
    elseif current_angle == 20
        Groups = [16;17;18];
    else
        Groups = find(Array_angles == current_angle);
    end

    for ii = 1: length(Groups) % different exp. groups (might be different days)

        no_Group = Groups(ii);

        the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
        theTruestorePath = storePath{no_Group};
        theTruestorePath = strrep(theTruestorePath,'results','results_plus');
        thefiles = dir(fullfile(theTruestorePath,'*.mat'));

        for file_ind = 1:length(thefiles) % load the cases

            filename = thefiles(file_ind).name;

            if contains(filename, 'PlusInfo_')

                load(fullfile(thefiles(1).folder, thefiles(file_ind).name));
                ContourL_all = xy.arclen_spl(Good_case_frm);
                ContourL = VicFc_Get_ContourLength(ContourL_all) * mag; % unit: um
                if ContourL < 30
                    for frm_ind = 1:size(Good_case_frm,2) % different snapshots

                        xy_ind = Good_case_frm(frm_ind); % index of the 'good' cases
                        fiber_xy_in_lattice = [fiber_xy_in_lattice; xy.fiber_xy_in_lattice{1, xy_ind}];
                        curvature = [curvature; xy.curvature{1, xy_ind}/mag]; % /mag: convert curvature unit to 1/um

                    end
                end
            end
        end
    end

    fiber_xy_in_lattice(fiber_xy_in_lattice > 1) = 1;
    fiber_xy_in_lattice(fiber_xy_in_lattice < 0) = 0;  % to make sure no error in accumarray
    [~,~,~,ind_x,ind_y] = histcounts2(fiber_xy_in_lattice(:,1), ...
        fiber_xy_in_lattice(:,2), Xedges, Yedges);
    binned = accumarray([ind_y,ind_x],curvature,[numel(Xedges)-1 numel(Xedges)-1],@VicMean,0);

    [x_plot,y_plot] = meshgrid(0.5/plot_resolution:1/plot_resolution:1-0.5/plot_resolution, ...
        0.5/plot_resolution:1/plot_resolution:1-0.5/plot_resolution);
    figure('color', 'w'); set(gcf, 'Position', [100 100 600 600]);

    pcolor(x_plot,y_plot, binned);
    axis equal; axis off
    cmocean('deep'); caxis([0 1]); %colorbar;

    viscircles([0 0; 0 1; 1 0; 1 1], mean(radii)*ones(4, 1)/mean([Ctr2Ctr_x, Ctr2Ctr_y]),'Color','r');
    xlim([0 1]); ylim([0 1]);

    f = gcf;
    savefig(f,['F:\Processing & Results\' ...
        'Actin Filaments in Porous Media\Figures\Curvature Distribution in Lattice\' ...
        'actin_curvature_lattice_angle_',num2str(current_angle),'_L_lower_30.fig'])

    set(f,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
        'Actin Filaments in Porous Media\Figures\Curvature Distribution in Lattice\' ...
        'actin_curvature_lattice_angle_',num2str(current_angle),'_L_lower_30.eps']);

    clearvars fiber_xy_in_lattice curvature

end

% another 20 degree:
no_Group = 28;
fiber_xy_in_lattice = [];
curvature = [];

the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
theTruestorePath = storePath{no_Group};
theTruestorePath = strrep(theTruestorePath,'results','results_plus');
thefiles = dir(fullfile(theTruestorePath,'*.mat'));
for file_ind = 1:length(thefiles) % load the cases

    filename = thefiles(file_ind).name;

    if contains(filename, 'PlusInfo_')

        load(fullfile(thefiles(1).folder, thefiles(file_ind).name));

        ContourL_all = xy.arclen_spl(Good_case_frm);
        ContourL = VicFc_Get_ContourLength(ContourL_all) * mag; % unit: um
        if ContourL < 30
            for frm_ind = 1:size(Good_case_frm,2) % different snapshots

                xy_ind = Good_case_frm(frm_ind); % index of the 'good' cases
                fiber_xy_in_lattice = [fiber_xy_in_lattice; xy.fiber_xy_in_lattice{1, xy_ind}];
                curvature = [curvature; xy.curvature{1, xy_ind}/mag]; % /mag: convert curvature unit to 1/um

            end
        end
    end
end

fiber_xy_in_lattice(fiber_xy_in_lattice > 1) = 1;
fiber_xy_in_lattice(fiber_xy_in_lattice < 0) = 0;  % to make sure no error in accumarray
[~,~,~,ind_x,ind_y] = histcounts2(fiber_xy_in_lattice(:,1), ...
    fiber_xy_in_lattice(:,2), Xedges, Yedges);
binned = accumarray([ind_y,ind_x],curvature,[numel(Xedges)-1 numel(Xedges)-1],@VicMean,0);

[x_plot,y_plot] = meshgrid(0.5/plot_resolution:1/plot_resolution:1-0.5/plot_resolution, ...
    0.5/plot_resolution:1/plot_resolution:1-0.5/plot_resolution);
figure('color', 'w'); set(gcf, 'Position', [100 100 600 600]);

pcolor(x_plot,y_plot, binned);
axis equal; axis off
cmocean('deep'); caxis([0 1]); %colorbar;

viscircles([0 0; 0 1; 1 0; 1 1], mean(radii)*ones(4, 1)/mean([Ctr2Ctr_x, Ctr2Ctr_y]),'Color','r');
xlim([0 1]); ylim([0 1]);

f = gcf;
savefig(f,['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Curvature Distribution in Lattice\' ...
    'actin_curvature_lattice_angle_20_smallerpillar_L_lower_30.fig'])

set(f,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Curvature Distribution in Lattice\' ...
    'actin_curvature_lattice_angle_20_smallerpillar_L_lower_30.eps']);



%% Plot fiber curvatures in the lattice (L < 30um; with colorbar).

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
Array_angles = cell2mat(xlsfile(:, 14));  % The flow angles.

plot_resolution = 50;
Xedges = [0:1/plot_resolution:1]; Yedges = [0:1/plot_resolution:1]; % for the histogram

for current_angle = [0 30 45] % different angles

    fiber_xy_in_lattice = [];
    curvature = [];

    if current_angle == 0
        Groups = [7;8];
    elseif current_angle == 20
        Groups = [16;17;18];
    else
        Groups = find(Array_angles == current_angle);
    end

    for ii = 1: length(Groups) % different exp. groups (might be different days)

        no_Group = Groups(ii);

        the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
        theTruestorePath = storePath{no_Group};
        theTruestorePath = strrep(theTruestorePath,'results','results_plus');
        thefiles = dir(fullfile(theTruestorePath,'*.mat'));

        for file_ind = 1:length(thefiles) % load the cases

            filename = thefiles(file_ind).name;

            if contains(filename, 'PlusInfo_')

                load(fullfile(thefiles(1).folder, thefiles(file_ind).name));
                ContourL_all = xy.arclen_spl(Good_case_frm);
                ContourL = VicFc_Get_ContourLength(ContourL_all) * mag; % unit: um
                if ContourL < 30
                    for frm_ind = 1:size(Good_case_frm,2) % different snapshots

                        xy_ind = Good_case_frm(frm_ind); % index of the 'good' cases
                        fiber_xy_in_lattice = [fiber_xy_in_lattice; xy.fiber_xy_in_lattice{1, xy_ind}];
                        curvature = [curvature; xy.curvature{1, xy_ind}/mag]; % /mag: convert curvature unit to 1/um

                    end
                end
            end
        end
    end

    fiber_xy_in_lattice(fiber_xy_in_lattice > 1) = 1;
    fiber_xy_in_lattice(fiber_xy_in_lattice < 0) = 0;  % to make sure no error in accumarray
    [~,~,~,ind_x,ind_y] = histcounts2(fiber_xy_in_lattice(:,1), ...
        fiber_xy_in_lattice(:,2), Xedges, Yedges);
    binned = accumarray([ind_y,ind_x],curvature,[numel(Xedges)-1 numel(Xedges)-1],@VicMean,0);

    [x_plot,y_plot] = meshgrid(0.5/plot_resolution:1/plot_resolution:1-0.5/plot_resolution, ...
        0.5/plot_resolution:1/plot_resolution:1-0.5/plot_resolution);
    figure('color', 'w'); set(gcf, 'Position', [100 100 800 600]);

    pcolor(x_plot,y_plot, binned);
    axis equal; axis off
    cmocean('deep'); caxis([0 1]); %colorbar;

    c = colorbar('Ticks',[0:0.5:1], 'Position',[0.83,0.20,0.03,0.60], ...
        'FontName','Times New Roman', 'FontSize', 24);  
    c.Label.String = '$\rm{Curvature}\, (1/\mu m)$';
    c.Label.Interpreter = 'LaTeX';
    c.TickLabelInterpreter = 'LaTeX';
    
    viscircles([0 0; 0 1; 1 0; 1 1], mean(radii)*ones(4, 1)/mean([Ctr2Ctr_x, Ctr2Ctr_y]),'Color','r');
    xlim([0 1]); ylim([0 1]);

    f = gcf;
    set(f,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
        'Actin Filaments in Porous Media\Figures\Curvature Distribution in Lattice\' ...
        'actin_curvature_lattice_angle_',num2str(current_angle),'_L_lower_30_colorbar.eps']);

    clearvars fiber_xy_in_lattice curvature

end


function y = VicMean(x)
xx = sort(x);
trash_no = floor(length(x) * 0.5);
y = mean(xx(trash_no+1: end));
end