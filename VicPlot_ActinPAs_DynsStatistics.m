%% Statistics about different actin dynamics in PAs.

clear; close all; clc;
xlsFolder = dir(['F:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Dynamics manually classification\*.xlsx']);

case_names = [];             ContourLs = [];             mu_bars = [];
FlowAngle = [];          
if_C_shape = [];        if_U_shape = [];            if_S_shape = [];
if_plus_S_shape = [];   if_buckled_3D = [];
if_folded = [];         if_coiled = [];             if_tumble = [];
if_has_note = [];       if_out_of_plane = [];       if_bad_reconstruction = [];

for ii = 1: length(xlsFolder)

    xlsfile = readcell([xlsFolder(1).folder, filesep, xlsFolder(ii).name], ...
        'Sheet','Sheet1','NumHeaderLines',1);

    case_names = [case_names; xlsfile(:, 1)];
    ContourLs = [ContourLs, xlsfile{:, 2}];
    mu_bars = [mu_bars, xlsfile{:, 14}];

    FlowAngle = [FlowAngle, xlsfile{:, 15}];

    if_C_shape = [if_C_shape, xlsfile{:, 3}];
    if_U_shape = [if_U_shape, xlsfile{:, 4}];
    if_S_shape = [if_S_shape, xlsfile{:, 5}];
    if_plus_S_shape = [if_plus_S_shape, xlsfile{:, 6}];

    if_buckled_3D = [if_buckled_3D, xlsfile{:, 7}];
    if_folded = [if_folded, xlsfile{:, 8}];
    if_coiled = [if_coiled, xlsfile{:, 9}];

    if_tumble = [if_tumble, xlsfile{:, 10}];

    if_has_note = [if_has_note, xlsfile{:, 11}];
    if_out_of_plane = [if_out_of_plane, xlsfile{:, 12}];
    if_bad_reconstruction = [if_bad_reconstruction, xlsfile{:, 13}];

end

f = figuresetting('centimeters',20,16,'times new roman',20,'on',1,'off','off');
f.figure('');
% for legends
semilogy(nan, nan, 'o','MarkerSize', 10,'MarkerEdgeColor','m','LineWidth',2); hold on
semilogy(nan, nan, '^','MarkerSize', 9,'MarkerEdgeColor','g','LineWidth',2); hold on
% plot
semilogy(FlowAngle(logical(if_C_shape)), mu_bars(logical(if_C_shape)), 'o','MarkerSize', 10,'MarkerEdgeColor','m','LineWidth',2); hold on
semilogy(FlowAngle(logical(if_U_shape)), mu_bars(logical(if_U_shape)), 'o','MarkerSize', 10,'MarkerEdgeColor','m','LineWidth',2); hold on
semilogy(FlowAngle(logical(if_S_shape)), mu_bars(logical(if_S_shape)), 'o','MarkerSize', 10,'MarkerEdgeColor','m','LineWidth',2); hold on
semilogy(FlowAngle(logical(if_plus_S_shape)), mu_bars(logical(if_plus_S_shape)), 'o','MarkerSize', 10,'MarkerEdgeColor','m','LineWidth',2); hold on
semilogy(FlowAngle(logical(if_buckled_3D)), mu_bars(logical(if_buckled_3D)), '^','MarkerSize', 9,'MarkerEdgeColor','g','LineWidth',2); hold on
semilogy(FlowAngle(logical(if_folded)), mu_bars(logical(if_folded)), '^','MarkerSize', 9,'MarkerEdgeColor','g','LineWidth',2); hold on
semilogy(FlowAngle(logical(if_coiled)), mu_bars(logical(if_coiled)), '^','MarkerSize', 9,'MarkerEdgeColor','g','LineWidth',2); hold on

f.interp_font('latex')
f.axes('linear',[-5 50],'log',[0 1e7],'$\alpha\,(^\circ)$','$\bar{\mu}$',24);
f.axes_ticks([0:5:45], 10.^[1:7]); grid on
legend({'2D deformation','3D deformation'}, 'FontName','Times New Roman', ...
    'Box','off','Location','southeast'); 

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Dynamics\actin_dyn_mubar.eps']);


% calculate mu_bar * flow_strength
FlowStrength = FlowAngle;
FlowStrength(FlowStrength == 0) = 3.2503;
FlowStrength(FlowStrength == 10) = 3.4097;
FlowStrength(FlowStrength == 15) = 3.2204;
FlowStrength(FlowStrength == 20) = 3.1802;
FlowStrength(FlowStrength == 30) = 2.8778;
FlowStrength(FlowStrength == 35) = 2.6439;
FlowStrength(FlowStrength == 45) = 2.3027;

mu_bars_times_flow_strength = mu_bars .* FlowStrength;

f = figuresetting('centimeters',20,16,'times new roman',20,'on',1,'off','off');
f.figure('');
% for legends
semilogy(nan, nan, 'o','MarkerSize', 10,'MarkerEdgeColor','m','LineWidth',2); hold on
semilogy(nan, nan, '^','MarkerSize', 9,'MarkerEdgeColor','g','LineWidth',2); hold on
% plot
semilogy(FlowAngle(logical(if_C_shape)), mu_bars_times_flow_strength(logical(if_C_shape)), ...
    'o','MarkerSize', 10,'MarkerEdgeColor','m','LineWidth',2); hold on
semilogy(FlowAngle(logical(if_U_shape)), mu_bars_times_flow_strength(logical(if_U_shape)), ...
    'o','MarkerSize', 10,'MarkerEdgeColor','m','LineWidth',2); hold on
semilogy(FlowAngle(logical(if_S_shape)), mu_bars_times_flow_strength(logical(if_S_shape)), ...
    'o','MarkerSize', 10,'MarkerEdgeColor','m','LineWidth',2); hold on
semilogy(FlowAngle(logical(if_plus_S_shape)), mu_bars_times_flow_strength(logical(if_plus_S_shape)), ...
    'o','MarkerSize', 10,'MarkerEdgeColor','m','LineWidth',2); hold on
semilogy(FlowAngle(logical(if_buckled_3D)), mu_bars_times_flow_strength(logical(if_buckled_3D)), ...
    '^','MarkerSize', 9,'MarkerEdgeColor','g','LineWidth',2); hold on
semilogy(FlowAngle(logical(if_folded)), mu_bars_times_flow_strength(logical(if_folded)), ...
    '^','MarkerSize', 9,'MarkerEdgeColor','g','LineWidth',2); hold on
semilogy(FlowAngle(logical(if_coiled)), mu_bars_times_flow_strength(logical(if_coiled)), ...
    '^','MarkerSize', 9,'MarkerEdgeColor','g','LineWidth',2); hold on

f.interp_font('latex')
f.axes('linear',[-5 50],'log',[0 5e7],'$\alpha\,(^\circ)$','$\bar{\mu} \delta$',24);
f.axes_ticks([0:5:45], 10.^[1:7]); grid on
legend({'2D deformation','3D deformation'}, 'FontName','Times New Roman', ...
    'Box','off','Location','southeast'); 

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Dynamics\actin_dyn_mubar_times_delta.eps']);


%% Statistics about actin-pillar dynamics in PAs.

clear; close all; clc;
xlsFolder = dir(['F:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Dynamics manually classification\*.xlsx']);

case_names = [];             ContourLs = [];             mu_bars = [];
FlowAngle = [];          
if_trapping = [];        if_passing = [];            if_others = [];

for ii = 1: length(xlsFolder)

    xlsfile = readcell([xlsFolder(1).folder, filesep, xlsFolder(ii).name], ...
        'Sheet','Sheet1','NumHeaderLines',1);

    case_names = [case_names; xlsfile(:, 1)];
    ContourLs = [ContourLs, xlsfile{:, 2}];
    mu_bars = [mu_bars, xlsfile{:, 14}];

    FlowAngle = [FlowAngle, xlsfile{:, 15}];

    if_trapping = [if_trapping, xlsfile{:, 16}];
    if_passing = [if_passing, xlsfile{:, 17}];
    if_others = [if_others, xlsfile{:, 18}];
 
end

FlowAngle_trapping = FlowAngle(logical(if_trapping));
FlowAngle_passing = FlowAngle(logical(if_passing));
FlowAngle_others = FlowAngle(logical(if_others));

ContourLs_trapping = ContourLs(logical(if_trapping));
ContourLs_passing = ContourLs(logical(if_passing));
ContourLs_others = ContourLs(logical(if_others));

f = figuresetting('centimeters',20,16,'times new roman',20,'on',1,'off','off');
f.figure('');
plot(FlowAngle_passing, ContourLs_passing, 'o','MarkerSize', 10,'MarkerEdgeColor','red','LineWidth',2); hold on
plot(FlowAngle_others, ContourLs_others, '^','MarkerSize', 9,'MarkerEdgeColor','blue','LineWidth',2); hold on
plot(FlowAngle_trapping, ContourLs_trapping, 'diamond','MarkerSize', 9,'MarkerEdgeColor',[0.92, 0.70, 0.22],'LineWidth', 2); hold on

f.interp_font('latex')
f.axes('linear',[-5 50],'linear',[0 80],'$\alpha\,(^\circ)$','$L\,(\rm{\mu m})$',24);
f.axes_ticks([0:5:45], [0:10:80]); grid on

legend({'Passing','Gliding','Trapping'}, 'FontName','Times New Roman', 'Box','off'); 

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Dynamics\actin_pillar_interaction.eps']);
