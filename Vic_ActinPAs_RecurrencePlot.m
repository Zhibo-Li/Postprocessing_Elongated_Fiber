% Poincaré plot or recurrence plot (RP).
% both experiment and simulation (tracers)


%% Calculte and save (actin filament in porous media)
clear; close all; clc;

xlsfile = readcell('ForActinPostprocessing.xlsx','Sheet','Sheet1','NumHeaderLines',1);
% This is the file that contains all the information about the later processing (in sheet 1).

NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.
PAsPath = xlsfile(:, 3);  % Path of the pillar array information.
Array_angles = xlsfile(:, 14);  % The flow angles.

n = 1;
for no_Group = [7 8 13:28]
    % Square-based array 0°, 10°, 20°, 15°, 35°, 30° and 45°

    Array_angle = Array_angles{no_Group};

    RotMatrix_init = rotz(-Array_angle); RotMatrix_init = RotMatrix_init(1:2, 1:2);
    % to rotate the pillar array

    the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
    filelist = dir(fullfile(storePath{no_Group},'*.mat'));
    % list of the .mat files which contain the reconstruction information
    % (came from 'Filaments detection' code) in one group.

    for no_Case = 1:length(filelist)

        filename = filelist(no_Case).name;
        if contains(filename, 'PAsInfoAdded_')
            load([storePath{no_Group}, filesep , filename])
            centers(:, 2) = 2048 - centers(:, 2);  % flip to image coordinate
            centers_new = (RotMatrix_init * centers')'; % rotate based on the design
            % viscircles(centers_new, radii,'LineStyle','--', 'LineWidth', 0.5, 'Color', 'k'); axis equal

            % correct the rotation degree
            sorted_centers_new = sortrows(centers_new,2);
            sorty = sort(centers_new(:, 2));
            sorty_ind = find(diff(sorty) > 100);
            corrected_angle = zeros(length(sorty_ind)-1, 1);
            for ii = 1: length(sorty_ind)-1
                to_be_fitted = sorted_centers_new(sorty_ind(ii)+1:sorty_ind(ii+1), :);
                fit_linear = fit(to_be_fitted(:, 1), to_be_fitted(:, 2), 'poly1');
                k = fit_linear.p1;
                corrected_angle(ii) = atand(k);
            end
            corrected_angle = mean(corrected_angle);
            RotMatrix_correct = rotz(-Array_angle-corrected_angle);
            RotMatrix_correct = RotMatrix_correct(1:2, 1:2);
            centers_corrected = (RotMatrix_correct * centers')';
            % viscircles(centers_corrected, radii,'LineStyle','--', 'LineWidth', 0.5, 'Color', 'r'); axis equal

            % get the approximate (after-gridlization) pillar center positions
            sortxx = sort(centers_corrected(:, 1)); sortyy = sort(centers_corrected(:, 2));
            sortxx_ind = find(diff(sortxx) > 100); sortyy_ind = find(diff(sortyy) > 100);
            PAs_X = zeros(length(sortxx_ind)+1, 1); PAs_Y = zeros(length(sortyy_ind)+1, 1);
            for jj = 1: length(sortxx_ind)-1
                PAs_X(jj+1) = mean(sortxx(sortxx_ind(jj)+1:sortxx_ind(jj+1)));
            end
            for jj = 1: length(sortyy_ind)-1
                PAs_Y(jj+1) = mean(sortyy(sortyy_ind(jj)+1:sortyy_ind(jj+1)));
            end
            PAs_X(1) = PAs_X(2) - mean(diff(PAs_X(2:end-1)));
            PAs_X(end) = PAs_X(end-1) + mean(diff(PAs_X(2:end-1)));
            PAs_Y(1) = PAs_Y(2) - mean(diff(PAs_Y(2:end-1)));
            PAs_Y(end) = PAs_Y(end-1) + mean(diff(PAs_Y(2:end-1)));
            % [XX, YY] = meshgrid(PAs_X, PAs_Y); plot(XX(:), YY(:))

            % average gap along y-direction
            ave_y_gap = mean(diff(PAs_Y));

            % get the filament length
            spl_Ls = xy.arclen_spl(Good_case_frm);
            L_0 = VicFc_Get_ContourLength(spl_Ls);

            % get the filament CoMs at different time
            fiber_center = reshape(cell2mat(xy.centroid), 2, length(cell2mat(xy.centroid))/2)';
            fiber_center = fiber_center(Good_case_frm, :);
            fiber_center = fiber_center + [lzero, lzero];
            fiber_center = (RotMatrix_correct * fiber_center')';
            % plot(fiber_center(:, 1), fiber_center(:, 2)); hold on
            % viscircles(centers_corrected, radii,'LineStyle','--', 'LineWidth', 0.5, 'Color', 'r'); axis equal

            % divide the trajectory into pieces according to their y-position
            fiber_Y_indicator = fiber_center(:, 2);
            for kk = 1: length(PAs_Y)-1
                fiber_Y_indicator(fiber_Y_indicator > PAs_Y(kk) & fiber_Y_indicator < PAs_Y(kk+1)) = kk;
            end

            % the 'entering lattice' positions
            ind_ToBeMoved = min(abs(repmat(fiber_center(:, 1), 1, length(PAs_X)) - PAs_X'), [], 2) > 60;
            fiber_center(ind_ToBeMoved, :) = [];  % only keep the cases that close to the lattice verticle edge
            fiber_Y_indicator(ind_ToBeMoved) = [];  % Y indicator as well
            fiber_X = fiber_center(:, 1);
            For_Poincare = nan(length(PAs_X), 3);
            for kk = 1: length(PAs_X)
                to_be_fitted2 = fiber_center(abs(fiber_X-PAs_X(kk))<=60, :);
                if ~isempty(to_be_fitted2) && numel(to_be_fitted2) > 2
                    fit_linear2 = fit(to_be_fitted2(:, 1), to_be_fitted2(:, 2), 'poly1');
                    For_Poincare(kk, 1) = mod((fit_linear2(PAs_X(kk))-PAs_Y(1)), ave_y_gap) / ave_y_gap;
                    % Correct: add '-PAs_Y(1)' @ 20230217 
                    For_Poincare(kk, 2) = kk;
                    For_Poincare(kk, 3) = fiber_Y_indicator(to_be_fitted2(1,1)==fiber_center(:,1));
                    % Correct: change 'fiber_Y_indicator(kk)' to fiber_Y_indicator(to_be_fitted2(1,1)==XXYY_Simu_i(:,1))
                end
            end
            For_Poincare(:, 3) = For_Poincare(:, 3) - min(For_Poincare(:, 3)) + 1;

            Info.L(n) = L_0;
            Info.map{n} = For_Poincare;
            Info.date(n) = the_exp_date;
            Info.name{n} = filename(14:end-4);
            Info.FlowAngle(n) = Array_angle;
            n = n + 1;

            % f=gcf;
            % exportgraphics(f,['D:\Dropbox\GitHub\tmp', filesep, num2str(the_exp_date), filename(25:end-11), '.png'],'Resolution',100)

        end
    end
end

save(['F:\Processing & Results\Actin Filaments in Porous Media\Figures\' ...
    'Poincare plots\Poincare_Map_data_202306.mat'], 'Info')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% draw Poincare Map for each flow angle (colorcode filament contour length)
clear; close all; clc;

load('F:\Processing & Results\Actin Filaments in Porous Media\Figures\Poincare plots\Poincare_Map_data_202306.mat');

Obj_Mag = 0.1; % um/pixel

theFlAng = Info.FlowAngle;
[C, ia, ic] = unique(theFlAng,'stable');

for ii = 1: length(C)

    Lattice_in_all = []; Lattice_out_all = []; L_toPlot_all = [];
    current_deg = C(ii);

    % the contour lengths
    L_all = Info.L  * Obj_Mag;

    % in & out position in a unit cell
    current_deg_index = find(theFlAng == current_deg);
    for jj = 1:length(current_deg_index)
        Lattice_in = Info.map{1, current_deg_index(jj)};

        Lattice_out = [[Lattice_in(2:end, 1);nan], Lattice_in(:, 2:3)];
        out_in_diff = Lattice_out(:, 2) - Lattice_in(:, 2);

        Lattice_in = Lattice_in(out_in_diff==0)';
        Lattice_out = Lattice_out(out_in_diff==0)';

%         figure('color', 'w'); set(gcf, 'Position', [100 100 500 500]);
%         plot(Lattice_in, Lattice_out, '.', 'LineStyle', 'none','MarkerSize', 13);
%         xlim([0 1]); ylim([0 1]); ax=gca; ax.FontSize = 15;
%         Info.name{jj}

        Lattice_in_all = [Lattice_in_all, Lattice_in];
        Lattice_out_all = [Lattice_out_all, Lattice_out];
        L_toPlot_all = [L_toPlot_all, L_all(current_deg_index(jj)) * ones(1, numel(Lattice_in))];

    end

    f = figuresetting('centimeters',10,10,'times new roman',20,'off',2,'off','off');
    f.figure('');
    scatter(Lattice_in_all, Lattice_out_all, 20, L_toPlot_all, 'Filled', 'o','MarkerEdgeColor','none');
    cmocean('amp'); caxis([0 50]);
%     hcb=colorbar; 
%     set(hcb,'TickLabelInterpreter','latex','Fontsize',20);
%     title(hcb,'$L(\mu{m})$','FontSize', 20,'Interpreter', 'latex');
    title_txt = strcat('$\alpha=', num2str(C(ii)), '^{\circ}$');

    f.interp_font('latex')
    f.axes('linear',[0 1],'linear',[0 1],'$\eta_{i}$','$\eta_{i+1}$',20);
    f.axes_ticks([0:0.2:1], [0:0.2:1]);
    grid on
    text(0.05, 0.9, title_txt,'FontSize', 20, ...
        'Interpreter', 'latex','BackgroundColor',[.7 .7 .7])

    set(gcf,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
        'Actin Filaments in Porous Media\Figures\Poincare plots\', ...
        title_txt(2:end-9),'_colorcode-L.eps']);

end

%% Plot only the colorbar for the above figures.
figure('color', 'w'); set(gcf, 'Position', [100 100 200 380]);

cmocean('amp');
caxis([0 50])
c = colorbar;
c.Label.String = '$L\ (\mathrm{\mu m})$';
c.Label.Interpreter = 'LaTeX';
c.TickLabelInterpreter = 'LaTeX';
c.FontSize = 20;
c.Location = 'west';
axis off

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\Actin Filaments in Porous Media\' ...
        'Figures\Poincare plots\The_colorbar_L.eps']);

%% Draw Poincare for tracer from simulation
clear; close all; clc;

Simu_deg = 25;  % Change this to choose the PAs angle (flow angle)
Simu_data = readmatrix(['F:\Simulation\202208_differentFlowangles_relatedto_' ...
    '0811exp_45deg\',num2str(Simu_deg),'deg\Data\Streamline_forPoincare_moreLines.csv']);
XXYY_Simu = Simu_data(1:end, 13:14);  % (x, y) of the streamline
IntegrationTime = Simu_data(1:end, 4);  % use to separate the different streamlines

% [pks,locs] = findpeaks(-diff(IntegrationTime), 'MinPeakHeight', 50);
locs = find(IntegrationTime==0);
ave_ele = size(XXYY_Simu, 1) / (length(locs)-1);  % average element number of each trajectory

RotMatrix = rotz(-Simu_deg); RotMatrix = RotMatrix(1:2, 1:2);
XXYY_Simu = (RotMatrix * XXYY_Simu')';  % after rotation

PAs_X = (-30:3:30)*1e-4; PAs_Y = (-30:3:30)*1e-4;  % pillar positions
ave_y_gap = 3e-4;

f = figuresetting('centimeters',10,10,'times new roman',22,'off',2,'off','off');
f.figure('');
for streamline_i = round((length(locs)-150)/2):length(locs)-round((length(locs)-150)/2)  
    % choose streamlines not too close to the sidewalls (sometimes needs change)

    XXYY_Simu_i = XXYY_Simu(locs(streamline_i):locs(streamline_i+1)-1, :);
%     plot(XXYY_Simu_i(:, 1), XXYY_Simu_i(:, 2))

    XXYY_Simu_i(XXYY_Simu_i(:, 1) > 0.0015, :) = []; 
    XXYY_Simu_i(XXYY_Simu_i(:, 1) < -0.0015, :) = [];  % remove the points out of the pillar array (sometimes needs change)
    % divide the trajectory into pieces according to their y-position
    fiber_Y_indicator = XXYY_Simu_i(:, 2);
    for kk = 1: length(PAs_Y)-1
        fiber_Y_indicator(fiber_Y_indicator > PAs_Y(kk) & fiber_Y_indicator < PAs_Y(kk+1)) = kk;
    end

    % the 'entering lattice' positions
    ind_ToBeMoved = min(abs(repmat(XXYY_Simu_i(:, 1), 1, length(PAs_X)) - PAs_X), [], 2) > 6e-6;
    XXYY_Simu_i(ind_ToBeMoved, :) = [];  % only keep the cases that close to the lattice verticle edge
    fiber_Y_indicator(ind_ToBeMoved) = [];  % Y indicator as well
    fiber_X = XXYY_Simu_i(:, 1);
    For_Poincare = nan(length(PAs_X), 3);
    for kk = 1: length(PAs_X)
        to_be_fitted2 = XXYY_Simu_i(abs(fiber_X-PAs_X(kk))<=6e-6, :);
        if ~isempty(to_be_fitted2) && numel(to_be_fitted2) > 2
            fit_linear2 = fit(to_be_fitted2(:, 1), to_be_fitted2(:, 2), 'poly1');
            For_Poincare(kk, 1) = mod((fit_linear2(PAs_X(kk))-PAs_Y(1)), ave_y_gap) / ave_y_gap;
            % Correct: add '-PAs_Y(1)' @ 20230217
            For_Poincare(kk, 2) = kk;
            For_Poincare(kk, 3) = fiber_Y_indicator(to_be_fitted2(1,1)==XXYY_Simu_i(:,1));
            % Correct: change 'fiber_Y_indicator(kk)' to fiber_Y_indicator(to_be_fitted2(1,1)==XXYY_Simu_i(:,1))
        end
    end

    For_Poincare(:, 3) = For_Poincare(:, 3) - min(For_Poincare(:, 3)) + 1;
    Lattice_in = For_Poincare;
    Lattice_out = [[Lattice_in(2:end, 1);nan], Lattice_in(:, 2:3)];
    out_in_diff = Lattice_out(:, 2) - Lattice_in(:, 2);

    Lattice_in = Lattice_in(out_in_diff==0)';
    Lattice_out = Lattice_out(out_in_diff==0)';

    plot(Lattice_in, Lattice_out, 'k.', 'LineStyle', 'none','MarkerSize', 8); hold on

end

f.interp_font('latex')
f.axes('linear',[0 1],'linear',[0 1],'$\eta_{i}=y_{\mathrm{in}}/\lambda$','$\eta_{i+1}=y_{\mathrm{out}}/\lambda$',22);
f.axes_ticks([0:0.2:1], [0:0.2:1]);
grid on
text(0.05, 0.9, ['$\alpha =',num2str(Simu_deg),'^{\circ}$'],'FontSize', 22, ...
    'Interpreter', 'latex','BackgroundColor',[.7 .7 .7])

f=gcf;
exportgraphics(f,['F:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Poincare plots\Simu_theta=',num2str(Simu_deg),'_tracer_newXYlabels.eps'])



%% Draw Poincare for tracer and actin (overlapping)
clear; close all; clc;

load('F:\Processing & Results\Actin Filaments in Porous Media\Figures\Poincare plots\Poincare_Map_data_202306.mat');

Obj_Mag = 0.1; % um/pixel

theFlAng = Info.FlowAngle;
[C, ia, ic] = unique(theFlAng,'stable');

for ii = 1: length(C)

    Lattice_in_all = []; Lattice_out_all = []; L_toPlot_all = [];
    current_deg = C(ii);

    % the contour lengths
    L_all = Info.L  * Obj_Mag;

    % in & out position in a unit cell
    current_deg_index = find(theFlAng == current_deg);
    for jj = 1:length(current_deg_index)
        Lattice_in = Info.map{1, current_deg_index(jj)};

        Lattice_out = [[Lattice_in(2:end, 1);nan], Lattice_in(:, 2:3)];
        out_in_diff = Lattice_out(:, 2) - Lattice_in(:, 2);

        Lattice_in = Lattice_in(out_in_diff==0)';
        Lattice_out = Lattice_out(out_in_diff==0)';

        Lattice_in_all = [Lattice_in_all, Lattice_in];
        Lattice_out_all = [Lattice_out_all, Lattice_out];
        L_toPlot_all = [L_toPlot_all, L_all(current_deg_index(jj)) * ones(1, numel(Lattice_in))];

    end

    f = figuresetting('centimeters',10,10,'times new roman',20,'off',2,'off','off');
    f.figure('');
    scatter(Lattice_in_all, Lattice_out_all, 20, L_toPlot_all, 'Filled', 'o','MarkerEdgeColor','none');
    cmocean('amp'); caxis([0 50]);

    hold on
    % plot for tracer
    Simu_deg = current_deg;  
    Simu_data = readmatrix(['F:\Simulation\202208_differentFlowangles_relatedto_' ...
        '0811exp_45deg\',num2str(Simu_deg),'deg\Data\Streamline_forPoincare_moreLines.csv']);
    XXYY_Simu = Simu_data(1:end, 13:14);  % (x, y) of the streamline
    IntegrationTime = Simu_data(1:end, 4);  % use to separate the different streamlines

    locs = find(IntegrationTime==0);
    ave_ele = size(XXYY_Simu, 1) / (length(locs)-1);  % average element number of each trajectory

    RotMatrix = rotz(-Simu_deg); RotMatrix = RotMatrix(1:2, 1:2);
    XXYY_Simu = (RotMatrix * XXYY_Simu')';  % after rotation

    PAs_X = (-30:3:30)*1e-4; PAs_Y = (-30:3:30)*1e-4;  % pillar positions
    ave_y_gap = 3e-4;

    for streamline_i = round((length(locs)-150)/2):3:length(locs)-round((length(locs)-150)/2)

        XXYY_Simu_i = XXYY_Simu(locs(streamline_i):locs(streamline_i+1)-1, :);

        XXYY_Simu_i(XXYY_Simu_i(:, 1) > 0.0015, :) = [];
        XXYY_Simu_i(XXYY_Simu_i(:, 1) < -0.0015, :) = [];  
        fiber_Y_indicator = XXYY_Simu_i(:, 2);
        for kk = 1: length(PAs_Y)-1
            fiber_Y_indicator(fiber_Y_indicator > PAs_Y(kk) & fiber_Y_indicator < PAs_Y(kk+1)) = kk;
        end

        % the 'entering lattice' positions
        ind_ToBeMoved = min(abs(repmat(XXYY_Simu_i(:, 1), 1, length(PAs_X)) - PAs_X), [], 2) > 6e-6;
        XXYY_Simu_i(ind_ToBeMoved, :) = [];  % only keep the cases that close to the lattice verticle edge
        fiber_Y_indicator(ind_ToBeMoved) = [];  % Y indicator as well
        fiber_X = XXYY_Simu_i(:, 1);
        For_Poincare = nan(length(PAs_X), 3);
        for kk = 1: length(PAs_X)
            to_be_fitted2 = XXYY_Simu_i(abs(fiber_X-PAs_X(kk))<=6e-6, :);
            if ~isempty(to_be_fitted2) && numel(to_be_fitted2) > 2
                fit_linear2 = fit(to_be_fitted2(:, 1), to_be_fitted2(:, 2), 'poly1');
                For_Poincare(kk, 1) = mod((fit_linear2(PAs_X(kk))-PAs_Y(1)), ave_y_gap) / ave_y_gap;
                % Correct: add '-PAs_Y(1)' @ 20230217
                For_Poincare(kk, 2) = kk;
                For_Poincare(kk, 3) = fiber_Y_indicator(to_be_fitted2(1,1)==XXYY_Simu_i(:,1));
                % Correct: change 'fiber_Y_indicator(kk)' to fiber_Y_indicator(to_be_fitted2(1,1)==XXYY_Simu_i(:,1))
            end
        end

        For_Poincare(:, 3) = For_Poincare(:, 3) - min(For_Poincare(:, 3)) + 1;
        Lattice_in = For_Poincare;
        Lattice_out = [[Lattice_in(2:end, 1);nan], Lattice_in(:, 2:3)];
        out_in_diff = Lattice_out(:, 2) - Lattice_in(:, 2);

        Lattice_in = Lattice_in(out_in_diff==0)';
        Lattice_out = Lattice_out(out_in_diff==0)';

        plot(Lattice_in, Lattice_out, 'cx', 'LineStyle', 'none','MarkerSize', 3); hold on

    end

    title_txt = strcat('$\alpha=', num2str(C(ii)), '^{\circ}$');
    f.interp_font('latex')
    f.axes('linear',[0 1],'linear',[0 1],'$\eta_{i}$','$\eta_{i+1}$',20);
    f.axes_ticks([0:0.2:1], [0:0.2:1]);
    grid on
    text(0.05, 0.9, title_txt,'FontSize', 20, ...
        'Interpreter', 'latex','BackgroundColor',[.7 .7 .7])

    set(gcf,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
        'Actin Filaments in Porous Media\Figures\Poincare plots\', ...
        title_txt(2:end-9),'_colorcode-L_with_tracer.eps']);

    close

end


