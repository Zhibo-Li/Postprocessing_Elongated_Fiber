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
for no_Group = [7 8 13:25]
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
            L_0 = Vic_Get_ave_cutExtreme(spl_Ls, 0.2);

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
    'Poincare plots\Poincare_Map_data.mat'], 'Info')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% draw Poincare Map for each flow angle (colorcode filament contour length)
clear; close all; clc;

load('F:\Processing & Results\Actin Filaments in Porous Media\Figures\Poincare plots\Poincare_Map_data.mat');

Obj_Mag = 0.1; % um/pixel

theFlAng = Info.FlowAngle;
[C, ia, ic] = unique(theFlAng,'stable');
ia = [ia; length(ic)+1];

L_all = [];
for ii = 1:length(ia)-1

    L_all = [L_all, Info.L(ia(ii):ia(ii+1)-1) * Obj_Mag];
    Lattice_in_all = []; Lattice_out_all = []; L_toPlot_all = [];

    proc_date = C(ii);
    for jj = ia(ii):ia(ii+1)-1
        Lattice_in = Info.map{1, jj};
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
        L_toPlot_all = [L_toPlot_all, L_all(jj) * ones(1, numel(Lattice_in))];

    end

    figure('color', 'w'); set(gcf, 'Position', [100 100 600 600]);
    scatter(Lattice_in_all, Lattice_out_all, 20, L_toPlot_all, 'Filled', 'o','MarkerEdgeColor','none');
    hcb=colorbar; caxis([0 50]); colormap jet
    set(hcb,'TickLabelInterpreter','latex','Fontsize',16);
    title(hcb,'$L(\mu{m})$','FontSize', 20,'Interpreter', 'latex');

    axis equal; grid on; box on
    xlim([0 1]); ylim([0 1]); 
    set(gca,'TickLabelInterpreter','latex','Fontsize',16);
    xlabel('$\eta_{i}$','FontSize', 24,'Interpreter', 'latex');
    ylabel('$\eta_{i+1}$','FontSize', 24,'Interpreter', 'latex');
    title_txt = strcat('$\theta=', num2str(C(ii)), '^{\circ}$');
    title(title_txt,'FontSize', 20,'Interpreter', 'latex');
    
%     f=gcf;
%     exportgraphics(f,['F:\Processing & Results\Actin Filaments in Porous Media\Figures\Poincare plots\',title_txt(2:end-9),'_colorcode-L.png'],'Resolution',100)

end

%% draw Poincare Map based on contour length (shorter and longer) at different flow angles
clear; close all; clc;

load('F:\Processing & Results\Actin Filaments in Porous Media\Figures\Poincare plots\Poincare_Map_data.mat');

Obj_Mag = 0.1; % um/pixel

theFlAng = Info.FlowAngle;
[C, ia, ic] = unique(theFlAng,'stable');
ia = [ia; length(ic)+1];

ii = 2; % choose theta = 10 
PAs_deg = 10;

proc_date = C(ii);
L_all = Info.L(ia(ii):ia(ii+1)-1);
Contour_L_divide = 150;
longer_ind = find(L_all >= Contour_L_divide);
shorter_ind = find(L_all < Contour_L_divide);

map_data = Info.map(1, ia(ii):ia(ii+1)-1);

%%%%% draw longer fiber 
Lattice_in_all = []; Lattice_out_all = [];
for jj = 1:length(longer_ind)

    Lattice_in = map_data{longer_ind(jj)};
    Lattice_out = [[Lattice_in(2:end, 1);nan], Lattice_in(:, 2:3)];
    out_in_diff = Lattice_out(:, 2) - Lattice_in(:, 2);

    Lattice_in = Lattice_in(out_in_diff==0)';
    Lattice_out = Lattice_out(out_in_diff==0)';

    Lattice_in_all = [Lattice_in_all, Lattice_in];
    Lattice_out_all = [Lattice_out_all, Lattice_out];
end

figure('color', 'w'); set(gcf, 'Position', [100 100 500 500]);
plot(Lattice_in_all, Lattice_out_all, '.', 'LineStyle', 'none','MarkerSize', 13);
axis equal; grid on
xlim([0 1]); ylim([0 1]); ax=gca; ax.FontSize = 15;
xlabel('$\eta_{i}$','FontSize', 22,'Interpreter', 'latex');
ylabel('$\eta_{i+1}$','FontSize', 22,'Interpreter', 'latex');
title(['$\theta=',num2str(PAs_deg),'^{\circ}:\ L>', num2str(Contour_L_divide * Obj_Mag),'\mu{m}$'],'FontSize', 20,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,['F:\Processing & Results\Actin Filaments in Porous Media\Figures\Poincare plots\theta=',num2str(PAs_deg),'_longer.png'],'Resolution',100)

%%%%% draw shorter fiber 
Lattice_in_all = []; Lattice_out_all = [];
for jj = 1:length(shorter_ind)

    Lattice_in = map_data{shorter_ind(jj)};
    Lattice_out = [[Lattice_in(2:end, 1);nan], Lattice_in(:, 2:3)];
    out_in_diff = Lattice_out(:, 2) - Lattice_in(:, 2);

    Lattice_in = Lattice_in(out_in_diff==0)';
    Lattice_out = Lattice_out(out_in_diff==0)';

    Lattice_in_all = [Lattice_in_all, Lattice_in];
    Lattice_out_all = [Lattice_out_all, Lattice_out];
end

figure('color', 'w'); set(gcf, 'Position', [100 100 500 500]);
plot(Lattice_in_all, Lattice_out_all, '.', 'LineStyle', 'none','MarkerSize', 13);
axis equal; grid on
xlim([0 1]); ylim([0 1]); ax=gca; ax.FontSize = 15;
xlabel('$\eta_{i}$','FontSize', 22,'Interpreter', 'latex');
ylabel('$\eta_{i+1}$','FontSize', 22,'Interpreter', 'latex');
title(['$\theta=',num2str(PAs_deg),'^{\circ}:\ L<', num2str(Contour_L_divide * Obj_Mag),'\mu{m}$'],'FontSize', 20,'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,['F:\Processing & Results\Actin Filaments in Porous Media\Figures\Poincare plots\theta=',num2str(PAs_deg),'_shorter.png'],'Resolution',100)



%% Draw Poincare for tracer from simulation
clear; close all; clc;

Simu_deg = 40;  % Change this to choose the PAs angle (flow angle)
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

figure('color', 'w'); set(gcf, 'Position', [100 100 500 500]);
for streamline_i = 250:length(locs)-249  % choose streamlines not too close to the sidewalls (sometimes needs change)

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

    plot(Lattice_in, Lattice_out, 'b.', 'LineStyle', 'none','MarkerSize', 13); hold on

end
ax=gca; ax.FontSize = 15; axis equal; grid on 
xlim([0 1]); ylim([0 1]); 
xlabel('$\eta_{i}$','FontSize', 22,'Interpreter', 'latex');
ylabel('$\eta_{i+1}$','FontSize', 22,'Interpreter', 'latex');
title(['$Simulation:\ \theta=',num2str(Simu_deg),'^{\circ}$'],'FontSize', 20,'Interpreter', 'latex');

f=gcf;
exportgraphics(f,['F:\Processing & Results\Actin Filaments in Porous Media\Figures\Poincare plots\Simu_theta=',num2str(Simu_deg),'_tracer.png'],'Resolution',100)


%% functions
function L_0 = Vic_Get_ave_cutExtreme(spl_Ls, threshold)
% the function is used to calculte the contour lenght of the filament based
% on the average without the extrame values.
%
% threshold -- threshold to truncate the extreme values

L_0_coarse = mean(spl_Ls, 'omitnan'); % roughly estimated filament length
spl_Ls(spl_Ls < L_0_coarse*(1-threshold)) = [];
spl_Ls(spl_Ls > L_0_coarse*(1+threshold)) = []; % remove the extreme value
if isempty(spl_Ls)
    L_0 = L_0_coarse;
    return
end
L_0 = mean(spl_Ls, 'omitnan');

while L_0_coarse ~= L_0 && ~isnan(L_0)
    L_0_coarse = L_0;
    spl_Ls(spl_Ls < L_0_coarse*(1-threshold)) = [];
    spl_Ls(spl_Ls > L_0_coarse*(1+threshold)) = []; % remove the extreme value
    if isempty(spl_Ls)
        L_0 = L_0_coarse;
        return
    end
    L_0 = mean(spl_Ls, 'omitnan');
end
end