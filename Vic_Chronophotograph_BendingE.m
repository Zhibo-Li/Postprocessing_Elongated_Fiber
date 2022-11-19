clear; close all; clc;

% % % % The .mat file where stores your results.
% [filename, pathname]=uigetfile({['F:\Processing & Results\Actin Filaments in Porous Media\' ...
%     '20210914-Actin\results\Group_1\*.mat']}, 'Choose a *.mat to be processed');  % input file
% load([pathname, filename])

% % % The .mat file where stores the PAs information.
load(['F:\Processing & Results\Actin Filaments in Porous Media\20210914-Actin' ...
    '\results\circlesforPAs20210914_G10.mat'])
centers(:, 2) = 2048 - centers(:, 2);

% % % The .mat files where stores your results.
thefiles = dir(['F:\Processing & Results\Actin Filaments in Porous Media\' ...
    'New calculations\20210914-Actin\results\Group_1\*.mat']);

for ii = 1: length(thefiles)

    load(fullfile(thefiles(1).folder, thefiles(ii).name));
    pathname = thefiles(1).folder;
    filename = thefiles(ii).name;

    figure('color', 'w'); set(gcf, 'Position', [100 300 1600 300]);
    Good_case_frm = find(ismember(xy(1).frame, Good_case));
%     Good_case_frm = Good_case;
    cmap = colormap('jet');
    for k = 1:size(Good_case_frm,2)

        j = Good_case_frm(k);% index of the 'good' cases
        lzero = 0;
%         try
%             lzero = max(lobject,ceil(5*lnoise)); % important!
%         catch
%             lzero = max(prmt(j).lobject,ceil(5*prmt(j).lnoise));
%         end

        if k == 1
            plot(xy.spl{j}(:,1)+lzero,xy.spl{j}(:,2)+lzero, 'LineWidth', 2);
            addaxislabel(1,'y (pixel)');
        elseif mod(k, 3) == 0
            addaxisplot(xy.spl{j}(:,1)+lzero,xy.spl{j}(:,2)+lzero, 1, 'color', ...
                cmap(mod(k*32, 255)+1, :), 'LineWidth', 2);
        end
        hold on

    end
    xlim([0 2050])
    xlabel('x (pixel)')
    axis equal; hold on

    B = 6.9e-26;  % Bending rigidity
    for k = 1:size(Good_case_frm,2)

        j = Good_case_frm(k);% index of the 'good' cases

        seglen_spl = xy.seglen_spl{1, j};
        % calculate the curvature
        X = xy(1).spl{j};
        [L2,R2,K2] = curvature(X);
        RR2 = movmean(R2,ceil(size(X,1)/100));

        XXX(k) = xy.centroid{1, j}(1);
        E(k) = B / 2 * sum((1./RR2(2:end)).^2 .* seglen_spl, 'omitnan') / 1e-7; % 1e-7 is the scale m/pixel;

    end

    addaxis(XXX+lzero, E, '*r', 'LineStyle','none', 'MarkerSize', 7);
    addaxislabel(2, 'Bending energy (J)');
    viscircles(centers, radii,'LineStyle','--', 'LineWidth', 0.5, 'Color', 'k'); hold on
    
    f=gcf;
    exportgraphics(f,[pathname, filesep, filename(1: end-4), '_BendingEnergy.png'],'Resolution',100)

    close all
    clearvars XXX E xy
end