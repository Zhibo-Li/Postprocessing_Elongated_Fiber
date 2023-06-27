%%% OpenFOAM data postprocessing: microchannel with the pillar array
% direction: x -- streamwise
%            y -- spanwise
%            z -- channel height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc

% Initial velocity
U_0 = 1e-5;

for Array_angle = 35

    % read simulation data
    data = readmatrix(['F:\Simulation\202208_differentFlowangles_relatedto_0811exp_45deg\', ...
        num2str(Array_angle), 'deg\Data\Deg', num2str(Array_angle), '_0o5Crop_Mid-Z0.csv']);
    Ux = data(1:end, 3);
    Uy = data(1:end, 4);
    % Uz = data(1:end, 5);
    XX = data(1:end, 6);
    YY = data(1:end, 7);
    % ZZ = data(1:end, 8);

    % Rotation matrix
    RotMatrix = rotz(-Array_angle); RotMatrix = RotMatrix(1:2, 1:2);

    % Rotate the velocity
    Rotated_velo = RotMatrix * [Ux, Uy]';
    Rotated_Ux = Rotated_velo(1, :);
    Rotated_Uy = Rotated_velo(2, :);

    % Rotate the position
    Rotated_pos = RotMatrix * [XX,YY]';
    Rotated_xx = Rotated_pos(1, :);
    Rotated_yy = Rotated_pos(2, :);

    % Interpolate to grid
    interpolant_Ux = scatteredInterpolant(Rotated_xx',Rotated_yy',Rotated_Ux','natural','none');
    interpolant_Uy = scatteredInterpolant(Rotated_xx',Rotated_yy',Rotated_Uy','natural','none');

    % Grid
    Grid_density = 1001;
    Grid_x1 = -6e-4; Grid_x2 = 6e-4; Grid_y1 = -6e-4; Grid_y2 = 8e-4;
    [xx,yy] = meshgrid(linspace(Grid_x1,Grid_x2,Grid_density), linspace(Grid_y1,Grid_y2,Grid_density));

    % Interpolate
    Ux_interp = interpolant_Ux(xx,yy);
    Uy_interp = interpolant_Uy(xx,yy);

    % To plot
    figure('color', 'w'); set(gcf, 'Position', [100 100 800 800]);
    % Plot contour: |U|
%     contourf(xx,yy,sqrt(Ux_interp.^2+Uy_interp.^2)/U_0,100,'LineStyle','none');
%     shading interp
    axis equal; axis off

    % Plot pillars
    hold on
    rectangle('Position',[-1 -1 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[-1 2 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[2 -1 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[2 2 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[-4 -1 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[-4 2 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[-4 -4 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[-1 -4 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');
    rectangle('Position',[2 -4 2 2]*1e-4, 'Curvature', [1 1], 'FaceColor', [.7 .7 .7], 'EdgeColor', 'none');

    xlim([-4 4]*1e-4)
    ylim([-4 4]*1e-4)

    ampcolor = cmocean('amp');
    densecolor = cmocean('dense');
%     caxis([0 5])

    % Plot streanlines
    hold on
    startX_ver_1 = -4.5e-4*ones(1, 50);
    startY_ver_1 = linspace(-6e-4,8e-4,50);
    lineobj = streamline(xx,yy,Ux_interp, Uy_interp, startX_ver_1(2:end-1), startY_ver_1(2:end-1));
    for ii = 1:numel(startX_ver_1)-2
        lineobj(ii).LineWidth = 1;
        lineobj(ii).Color = 'k';
        lineobj(ii).LineStyle = ':';
    end
    hold on
    
    % load the case 1
    load(['F:\Processing & Results\Actin Filaments in Porous Media\20220902-Actin' ...
        '\results_plus\Group_2\PlusInfo_trajectory_M63_Phi20_G10_FlAng35_0.5nM' ...
        '_1nL_Expo20ms_11_no19-interval2-no226_AABGR_batch1.mat'])

    % case 1 filament 1
    ind = ismember(xy(1).frame, 59);
    the_chose_fiber = xy.fiber_xy_in_lattice{1, ind};
    % x-direction (connect the filament)
    x_diff = diff(the_chose_fiber(:,1));
    if ~isempty(find(x_diff > 0.9, 1))
        cut_loc_x = find(x_diff > 0.9)'; cut_loc_x = [0, cut_loc_x, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_x) - 1
            the_chose_fiber(cut_loc_x(foo)+1: cut_loc_x(foo+1), 1) = ...
                the_chose_fiber(cut_loc_x(foo)+1: cut_loc_x(foo+1), 1) - 1*(foo - 1);
        end
    end
    if ~isempty(find(x_diff < -0.9, 1))
        cut_loc_x = find(x_diff  < -0.9)'; cut_loc_x = [0, cut_loc_x, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_x) - 1
            the_chose_fiber(cut_loc_x(foo)+1: cut_loc_x(foo+1), 1) = ...
                the_chose_fiber(cut_loc_x(foo)+1: cut_loc_x(foo+1), 1) + 1*(foo - 1);
        end
    end
    % y-direction (connect the filament)
    y_diff = diff(the_chose_fiber(:,2));
    if ~isempty(find(y_diff > 0.9, 1))
        cut_loc_y = find(y_diff > 0.9)'; cut_loc_y = [0, cut_loc_y, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_y) - 1
            the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) = ...
                the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) - 1*(foo - 1);
        end
    end
    if ~isempty(find(y_diff < -0.9, 1))
        cut_loc_y = find(y_diff < -0.9)'; cut_loc_y = [0, cut_loc_y, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_y) - 1
            the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) = ...
                the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) + 1*(foo - 1);
        end
    end
    plot((the_chose_fiber(:,1)) * 3.5e-4 - 2.9e-4, (the_chose_fiber(:,2) - 1) * 3.5e-4, ...
        'color', ampcolor(80, :), 'LineWidth', 5);  

    % case 1 filament 2
    ind = ismember(xy(1).frame, 67);
    the_chose_fiber = xy.fiber_xy_in_lattice{1, ind};
    % x-direction (connect the filament)
    x_diff = diff(the_chose_fiber(:,1));
    if ~isempty(find(x_diff > 0.9, 1))
        cut_loc_x = find(x_diff > 0.9)'; cut_loc_x = [0, cut_loc_x, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_x) - 1
            the_chose_fiber(cut_loc_x(foo)+1: cut_loc_x(foo+1), 1) = ...
                the_chose_fiber(cut_loc_x(foo)+1: cut_loc_x(foo+1), 1) - 1*(foo - 1);
        end
    end
    if ~isempty(find(x_diff < -0.9, 1))
        cut_loc_x = find(x_diff  < -0.9)'; cut_loc_x = [0, cut_loc_x, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_x) - 1
            the_chose_fiber(cut_loc_x(foo)+1: cut_loc_x(foo+1), 1) = ...
                the_chose_fiber(cut_loc_x(foo)+1: cut_loc_x(foo+1), 1) + 1*(foo - 1);
        end
    end
    % y-direction (connect the filament)
    y_diff = diff(the_chose_fiber(:,2));
    if ~isempty(find(y_diff > 0.9, 1))
        cut_loc_y = find(y_diff > 0.9)'; cut_loc_y = [0, cut_loc_y, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_y) - 1
            the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) = ...
                the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) - 1*(foo - 1);
        end
    end
    if ~isempty(find(y_diff < -0.9, 1))
        cut_loc_y = find(y_diff < -0.9)'; cut_loc_y = [0, cut_loc_y, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_y) - 1
            the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) = ...
                the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) + 1*(foo - 1);
        end
    end
    plot((the_chose_fiber(:,1) + 1) * 3.5e-4 - 2.9e-4, (the_chose_fiber(:,2) - 1) * 3.5e-4, ...
        'color', ampcolor(140, :), 'LineWidth', 5);  

    % case 1 filament 3
    ind = ismember(xy(1).frame, 78);
    the_chose_fiber = xy.fiber_xy_in_lattice{1, ind};
    % x-direction (connect the filament)
    x_diff = diff(the_chose_fiber(:,1));
    if ~isempty(find(x_diff > 0.9, 1))
        cut_loc_x = find(x_diff > 0.9)'; cut_loc_x = [0, cut_loc_x, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_x) - 1
            the_chose_fiber(cut_loc_x(foo)+1: cut_loc_x(foo+1), 1) = ...
                the_chose_fiber(cut_loc_x(foo)+1: cut_loc_x(foo+1), 1) - 1*(foo - 1);
        end
    end
    if ~isempty(find(x_diff < -0.9, 1))
        cut_loc_x = find(x_diff  < -0.9)'; cut_loc_x = [0, cut_loc_x, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_x) - 1
            the_chose_fiber(cut_loc_x(foo)+1: cut_loc_x(foo+1), 1) = ...
                the_chose_fiber(cut_loc_x(foo)+1: cut_loc_x(foo+1), 1) + 1*(foo - 1);
        end
    end
    % y-direction (connect the filament)
    y_diff = diff(the_chose_fiber(:,2));
    if ~isempty(find(y_diff > 0.9, 1))
        cut_loc_y = find(y_diff > 0.9)'; cut_loc_y = [0, cut_loc_y, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_y) - 1
            the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) = ...
                the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) - 1*(foo - 1);
        end
    end
    if ~isempty(find(y_diff < -0.9, 1))
        cut_loc_y = find(y_diff < -0.9)'; cut_loc_y = [0, cut_loc_y, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_y) - 1
            the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) = ...
                the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) + 1*(foo - 1);
        end
    end
    plot(the_chose_fiber(:,1) * 3.5e-4 + 3.1e-4, the_chose_fiber(:,2) * 3.5e-4 - 3e-4, ...
        'color', ampcolor(200, :), 'LineWidth', 5);  


    % load the case 2
    load(['F:\Processing & Results\Actin Filaments in Porous Media\20220902-Actin' ...
        '\results_plus\Group_1\PlusInfo_trajectory_M63_Phi20_G10_FlAng35_0.5nM' ...
        '_0.5nL_Expo20ms_31_no1-interval2-no291_AABGR_batch1.mat'])
    % case 2 filament 1
    ind = ismember(xy(1).frame, 49);
    the_chose_fiber = xy.fiber_xy_in_lattice{1, ind};
        % y-direction (connect the filament)
    y_diff = diff(the_chose_fiber(:,2));
    if ~isempty(find(y_diff > 0.9, 1))
        cut_loc_y = find(y_diff > 0.9)'; cut_loc_y = [0, cut_loc_y, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_y) - 1
            the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) = ...
                the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) - 1*(foo - 1);
        end
    end
    if ~isempty(find(y_diff < -0.9, 1))
        cut_loc_y = find(y_diff < -0.9)'; cut_loc_y = [0, cut_loc_y, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_y) - 1
            the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) = ...
                the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) + 1*(foo - 1);
        end
    end
    plot((the_chose_fiber(:,1)) * 3.5e-4 - 3e-4, the_chose_fiber(:,2) * 3.5e-4, ...
        'color', densecolor(80, :), 'LineWidth', 5);  

    % case 2 filament 2
    ind = ismember(xy(1).frame, 85);
    the_chose_fiber = xy.fiber_xy_in_lattice{1, ind};
        % y-direction (connect the filament)
    y_diff = diff(the_chose_fiber(:,2));
    if ~isempty(find(y_diff > 0.9, 1))
        cut_loc_y = find(y_diff > 0.9)'; cut_loc_y = [0, cut_loc_y, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_y) - 1
            the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) = ...
                the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) - 1*(foo - 1);
        end
    end
    if ~isempty(find(y_diff < -0.9, 1))
        cut_loc_y = find(y_diff < -0.9)'; cut_loc_y = [0, cut_loc_y, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_y) - 1
            the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) = ...
                the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) + 1*(foo - 1);
        end
    end
    plot((the_chose_fiber(:,1)) * 3.5e-4 - 3e-4, the_chose_fiber(:,2) * 3.5e-4 - 3e-4, ...
        'color', densecolor(140, :), 'LineWidth', 5); 

    % case 2 filament 3
    ind = ismember(xy(1).frame, 94);
    the_chose_fiber = xy.fiber_xy_in_lattice{1, ind};
        % y-direction (connect the filament)
    y_diff = diff(the_chose_fiber(:,2));
    if ~isempty(find(y_diff > 0.9, 1))
        cut_loc_y = find(y_diff > 0.9)'; cut_loc_y = [0, cut_loc_y, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_y) - 1
            the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) = ...
                the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) - 1*(foo - 1);
        end
    end
    if ~isempty(find(y_diff < -0.9, 1))
        cut_loc_y = find(y_diff < -0.9)'; cut_loc_y = [0, cut_loc_y, size(the_chose_fiber, 1)];
        for foo = 1: length(cut_loc_y) - 1
            the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) = ...
                the_chose_fiber(cut_loc_y(foo)+1: cut_loc_y(foo+1), 2) + 1*(foo - 1);
        end
    end
    plot((the_chose_fiber(:,1)) * 3.5e-4, the_chose_fiber(:,2) * 3.5e-4 - 3e-4, ...
        'color', densecolor(200, :), 'LineWidth', 5); 

    % draw a box
    plot([-4 4 4 -4 -4]*1e-4, [-4 -4 4 4 -4]*1e-4, 'k', 'LineWidth', 2);

    f=gcf;
    savefig(f,['D:\Dropbox\Research\My PhD thesis\Figures\5-flexible_fiber_array' ...
        '\DLD_fiber_sketch\DLD_fiber_sketch.fig'])
    set(f,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',['D:\Dropbox\Research\My PhD thesis' ...
        '\Figures\5-flexible_fiber_array\DLD_fiber_sketch\DLD_fiber_sketch.eps'])

end