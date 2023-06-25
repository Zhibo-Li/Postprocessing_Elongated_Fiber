%%%% Plot normalized the fiber x-y coordinates in the lattice.
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

Xedges = [0:0.02:1]; Yedges = [0:0.02:1]; % for the histogram

for current_angle = [0 10 15 20 30 35 45] % different angles

    fiber_xy_in_lattice = [];

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

                for frm_ind = 1:size(Good_case_frm,2) % different snapshots

                    xy_ind = Good_case_frm(frm_ind); % index of the 'good' cases
                    fiber_xy_in_lattice = [fiber_xy_in_lattice; xy.fiber_xy_in_lattice{1, xy_ind}];

                end
            end
        end
    end

    figure('color', 'w'); set(gcf, 'Position', [100 100 600 600]);
    h = histogram2(fiber_xy_in_lattice(:,1),fiber_xy_in_lattice(:,2), Xedges, ...
        Yedges, 'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
    axis equal; axis off
    cmocean('solar'); caxis([0 5]); %colorbar;
    
    viscircles([0 0; 0 1; 1 0; 1 1], mean(radii)*ones(4, 1)/mean([Ctr2Ctr_x, Ctr2Ctr_y]),'Color','r');
    xlim([0 1]); ylim([0 1]);

    f = gcf;
    savefig(f,['F:\Processing & Results\' ...
        'Actin Filaments in Porous Media\Figures\Position Distribution in Lattice\' ...
        'actin_pos_distri_lattice_angle_',num2str(current_angle),'.fig'])

    set(f,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
        'Actin Filaments in Porous Media\Figures\Position Distribution in Lattice\' ...
        'actin_pos_distri_lattice_angle_',num2str(current_angle),'.eps']);

    clearvars fiber_xy_in_lattice

end

% another 20 degree:
no_Group = 28;
fiber_xy_in_lattice = [];

the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
theTruestorePath = storePath{no_Group};
theTruestorePath = strrep(theTruestorePath,'results','results_plus');
thefiles = dir(fullfile(theTruestorePath,'*.mat'));
for file_ind = 1:length(thefiles) % load the cases

    filename = thefiles(file_ind).name;

    if contains(filename, 'PlusInfo_')

        load(fullfile(thefiles(1).folder, thefiles(file_ind).name));

        for frm_ind = 1:size(Good_case_frm,2) % different snapshots

            xy_ind = Good_case_frm(frm_ind); % index of the 'good' cases
            fiber_xy_in_lattice = [fiber_xy_in_lattice; xy.fiber_xy_in_lattice{1, xy_ind}];

        end
    end
end

figure('color', 'w'); set(gcf, 'Position', [100 100 600 600]);
h = histogram2(fiber_xy_in_lattice(:,1),fiber_xy_in_lattice(:,2), Xedges, ...
    Yedges, 'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
axis equal; axis off
cmocean('solar'); caxis([0 5]); %colorbar;

viscircles([0 0; 0 1; 1 0; 1 1], mean(radii)*ones(4, 1)/mean([Ctr2Ctr_x, Ctr2Ctr_y]),'Color','r');
xlim([0 1]); ylim([0 1]);

f = gcf;
savefig(f,['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Position Distribution in Lattice\' ...
    'actin_pos_distri_lattice_angle_20_smallerpillar.fig'])

set(f,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Position Distribution in Lattice\' ...
    'actin_pos_distri_lattice_angle_20_smallerpillar.eps']);



%% Plot only the colorbar for the above figures.
figure('color', 'w'); set(gcf, 'Position', [100 100 200 380]);

cmocean('solar');
caxis([0 5])
c = colorbar;
c.Label.String = '$\rm{PDF}$';
c.Label.Interpreter = 'LaTeX';
c.TickLabelInterpreter = 'LaTeX';
c.FontSize = 20;
c.Location = 'west';
axis off

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
        'Actin Filaments in Porous Media\Figures\Position Distribution in Lattice\' ...
        'actin_pos_distri_lattice_ColorBar.eps']);



%%
%%%% Plot normalized the fiber x-y coordinates in the lattice.
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

Xedges = [0:0.02:1]; Yedges = [0:0.02:1]; % for the histogram

for current_angle = [0 10 15 20 30 35 45] % different angles

    fiber_xy_in_lattice = [];

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

                    end
                end
            end
        end
    end

    figure('color', 'w'); set(gcf, 'Position', [100 100 600 600]);
    h = histogram2(fiber_xy_in_lattice(:,1),fiber_xy_in_lattice(:,2), Xedges, ...
        Yedges, 'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
    axis equal; axis off
    cmocean('solar'); caxis([0 5]); %colorbar;
    
    viscircles([0 0; 0 1; 1 0; 1 1], mean(radii)*ones(4, 1)/mean([Ctr2Ctr_x, Ctr2Ctr_y]),'Color','r');
    xlim([0 1]); ylim([0 1]);

    f = gcf;
    savefig(f,['F:\Processing & Results\' ...
        'Actin Filaments in Porous Media\Figures\Position Distribution in Lattice\' ...
        'actin_pos_distri_lattice_angle_',num2str(current_angle),'_L_lower_30.fig'])

    set(f,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
        'Actin Filaments in Porous Media\Figures\Position Distribution in Lattice\' ...
        'actin_pos_distri_lattice_angle_',num2str(current_angle),'_L_lower_30.eps']);

    clearvars fiber_xy_in_lattice

end

% another 20 degree:
no_Group = 28;
fiber_xy_in_lattice = [];

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

            end
        end
    end
end

figure('color', 'w'); set(gcf, 'Position', [100 100 600 600]);
h = histogram2(fiber_xy_in_lattice(:,1),fiber_xy_in_lattice(:,2), Xedges, ...
    Yedges, 'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
axis equal; axis off
cmocean('solar'); caxis([0 5]); %colorbar;

viscircles([0 0; 0 1; 1 0; 1 1], mean(radii)*ones(4, 1)/mean([Ctr2Ctr_x, Ctr2Ctr_y]),'Color','r');
xlim([0 1]); ylim([0 1]);

f = gcf;
savefig(f,['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Position Distribution in Lattice\' ...
    'actin_pos_distri_lattice_angle_20_smallerpillar_L_lower_30.fig'])

set(f,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Position Distribution in Lattice\' ...
    'actin_pos_distri_lattice_angle_20_smallerpillar_L_lower_30.eps']);



%%
%%%% Plot normalized the fiber x-y coordinates in the lattice.
% With colorbar

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

Xedges = [0:0.02:1]; Yedges = [0:0.02:1]; % for the histogram

for current_angle = [0 30 45] % different angles

    fiber_xy_in_lattice = [];

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

                    end
                end
            end
        end
    end

    figure('color', 'w'); set(gcf, 'Position', [100 100 700 600]);
    h = histogram2(fiber_xy_in_lattice(:,1),fiber_xy_in_lattice(:,2), Xedges, ...
        Yedges, 'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
    axis equal; axis off
    cmocean('solar'); 
    if current_angle == 0
        caxis([0 8]); colorbar('Ticks',[0:4:8], 'Position',[0.87,0.20,0.03,0.60], ...
            'FontName','Times New Roman', 'FontSize',32);  % angle = 0
    elseif current_angle == 30
        caxis([0 8]); colorbar('Ticks',[0:4:8], 'Position',[0.87,0.20,0.03,0.60], ...
            'FontName','Times New Roman', 'FontSize',32);  % angle = 30
    elseif current_angle == 45
        caxis([0 6]); colorbar('Ticks',[0:3:6], 'Position',[0.87,0.20,0.03,0.60], ...
            'FontName','Times New Roman', 'FontSize',32);  % angle = 45
    end
    
    viscircles([0 0; 0 1; 1 0; 1 1], mean(radii)*ones(4, 1)/mean([Ctr2Ctr_x, Ctr2Ctr_y]),'Color','r');
    xlim([0 1]); ylim([0 1]);

    f = gcf;
    set(f,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
        'Actin Filaments in Porous Media\Figures\Position Distribution in Lattice\' ...
        'actin_pos_distri_lattice_angle_',num2str(current_angle),'_L_lower_30_colorbar.eps']);

    clearvars fiber_xy_in_lattice

end



