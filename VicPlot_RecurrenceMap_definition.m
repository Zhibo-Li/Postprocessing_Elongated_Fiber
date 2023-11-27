%% Draw Poincare for tracer from simulation
clear; close all; clc;

Simu_deg = 20;  % Change this to choose the PAs angle (flow angle)
Simu_data = readmatrix(['D:\Dropbox\Transfer\202208_differentFlowangles_related' ...
    'to_0811exp_45deg\',num2str(Simu_deg),'deg\Data\Streamline_forPoincare_moreLines.csv']);
XXYY_Simu = Simu_data(1:end, 13:14);  % (x, y) of the streamline
IntegrationTime = Simu_data(1:end, 4);  % use to separate the different streamlines

% [pks,locs] = findpeaks(-diff(IntegrationTime), 'MinPeakHeight', 50);
locs = find(IntegrationTime==0);
ave_ele = size(XXYY_Simu, 1) / (length(locs)-1);  % average element number of each trajectory

RotMatrix = rotz(-Simu_deg); RotMatrix = RotMatrix(1:2, 1:2);
XXYY_Simu = (RotMatrix * XXYY_Simu')';  % after rotation

PAs_X = (-15:3:3)*1e-4; PAs_Y = (-12:3:-3)*1e-4;  % pillar positions
ave_y_gap = 3e-4;

streamline_i = 87;
% choose streamlines not too close to the sidewalls (sometimes needs change)

XXYY_Simu_i = XXYY_Simu(locs(streamline_i):locs(streamline_i+1)-1, :);
XXYY_Simu_i(XXYY_Simu_i(:, 1) > 0.0004, :) = [];
XXYY_Simu_i(XXYY_Simu_i(:, 1) < -0.0016, :) = [];  % remove the points out of the pillar array (sometimes needs change)

n = 1;
for x_cut = -0.0015:0.0003:0.0003
    inter_P = InterX(XXYY_Simu_i', [x_cut x_cut; -12e-4 -3e-4]);
    Lattice_in(n) = mod(inter_P(2), 0.0003)/0.0003;
    n = n + 1;
end

Lattice_out = [Lattice_in(2:end), nan];

for jj = 1:6

    seg_XXYY_Simu_i = XXYY_Simu_i;
    seg_XXYY_Simu_i(seg_XXYY_Simu_i(:, 1) > jj*0.0003-0.0015, :) = [];
    seg_XXYY_Simu_i(seg_XXYY_Simu_i(:, 1) < (jj-1)*0.0003-0.0015, :) = [];  % remove the points out of the pillar array (sometimes needs change)

    figure('color', 'w');
    set(gcf, 'Position',[100 100 1000 600], 'Color','white', 'DefaultTextInterpreter', 'latex')
    plot(XXYY_Simu_i(:, 1), XXYY_Simu_i(:, 2),'k--', 'linewidth', 1.5)
    hold on
    plot(seg_XXYY_Simu_i(:, 1), seg_XXYY_Simu_i(:, 2),'m', 'linewidth', 3)
    hold on
    [XXYY, YYXX] = meshgrid(PAs_X, PAs_Y);
    pillar_R = 1e-4; radii = pillar_R*ones(length(YYXX(:)), 1);
    viscircles([XXYY(:), YYXX(:)], radii, 'color', 'r', 'linewidth', 0.5);
    axis equal; axis off
    hold on
    rectangle('Position', [-0.0016 -0.0013 0.0020 0.0011],'LineWidth', 0.3, 'EdgeColor', [.1 .1 .1])

    set(gcf,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',['D:\Dropbox\Research\Zhibo PhD thesis' ...
        '\Defense\Figures and videos\Poincare\trajectory_',num2str(jj),'.eps']);

    f = figuresetting('centimeters',10,10,'times new roman',22,'off',2,'off','off');
    f.figure('');
    plot(Lattice_in(1:jj), Lattice_out(1:jj), 'm.', 'LineStyle', 'none','MarkerSize', 25); hold on

    f.interp_font('latex')
    f.axes('linear',[0 1],'linear',[0 1],'$\eta_{i}$','$\eta_{i+1}$',22);
    f.axes_ticks([0:0.2:1], [0:0.2:1]);
    grid on
    % text(0.05, 0.9, ['$\alpha =',num2str(Simu_deg),'^{\circ}$'],'FontSize', 22, ...
    %     'Interpreter', 'latex','BackgroundColor',[.7 .7 .7])

    set(gcf,'renderer','Painters');
    print('-depsc2','-tiff','-r100','-vector',['D:\Dropbox\Research\Zhibo PhD thesis' ...
        '\Defense\Figures and videos\Poincare\Poincare_',num2str(jj),'.eps']);

end



