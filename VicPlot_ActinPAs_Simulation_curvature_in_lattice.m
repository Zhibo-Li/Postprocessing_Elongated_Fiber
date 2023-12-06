%% Plot fiber curvatures in the lattice for simulation data.
% data from ..\simulations_results\Flow_Ang_alpha0_*\*.mat

clear; close all; clc;

set(0,'DefaultAxesFontSize',14);
set(0,'defaulttextfontsize',14);
% set(0,'defaultAxesFontName', 'times new roman');
% set(0,'defaultTextFontName', 'times new roman');
set(0,'defaultAxesFontName', 'Arial');
set(0,'defaultTextFontName', 'Arial');
set(0,'defaulttextInterpreter','latex')

plot_resolution = 50;
Xedges = [0:1/plot_resolution:1]; Yedges = [0:1/plot_resolution:1]; % for the histogram

parent_path = ['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\' ...
    'FSI - Actin in PAs\send_zhibo\simulations_results'];

for current_angle = 0:5:45 % different angles

    fiber_xy_in_lattice = [];
    curvature = [];

    current_load_path = [parent_path, filesep, 'Flow_Ang_alpha0_', num2str(current_angle), filesep];
    load_results_list = dir(fullfile(current_load_path, '\*.mat'));

    for ii = 1:length(load_results_list)
        load(fullfile(load_results_list(ii).folder, load_results_list(ii).name));

        for jj = 1:length(fiberInfo)
            fiber_xy_in_lattice = [fiber_xy_in_lattice; fiberInfo(jj).fiber_xy_in_lattice];
            curvature = [curvature; fiberInfo(jj).curvature]; % unit: 1/um
        end
    end

    fiber_xy_in_lattice(:,2) = 1 - fiber_xy_in_lattice(:,2); % flip vertically

    [~,~,~,ind_x,ind_y] = histcounts2(fiber_xy_in_lattice(:,1), ...
        fiber_xy_in_lattice(:,2), Xedges, Yedges);
    binned = accumarray([ind_y,ind_x],curvature,[numel(Xedges)-1 numel(Xedges)-1],@mean,0);

    [x_plot,y_plot] = meshgrid(0.5/plot_resolution:1/plot_resolution:1-0.5/plot_resolution, ...
        0.5/plot_resolution:1/plot_resolution:1-0.5/plot_resolution);
    figure('color', 'w'); set(gcf, 'Position', [100 100 600 600]);
    
    pcolor(x_plot,y_plot, binned);
    axis equal; axis off
    cmocean('deep'); caxis([0 1]); %colorbar;
    
    viscircles([0 0; 0 1; 1 0; 1 1], 1/3*ones(4, 1),'Color','r');
    xlim([0 1]); ylim([0 1]);

    f = gcf;
    savefig(f,['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\FSI - Actin in PAs' ...
        '\send_zhibo\figures\Curvature Distribution in Lattice\' ...
        'actin_curvature_lattice_angle_',num2str(current_angle),'_Simulation.fig'])

    set(f,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',['D:\Dropbox\Collaboration - LadHyX\' ...
        'Give_to_Zhibo_nonShared\FSI - Actin in PAs\send_zhibo\figures\Curvature Distribution in Lattice\' ...
        'actin_curvature_lattice_angle_',num2str(current_angle),'_Simulation.eps']);

    clearvars fiber_xy_in_lattice curvature

end