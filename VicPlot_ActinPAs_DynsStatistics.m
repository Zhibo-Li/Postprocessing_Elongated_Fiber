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
if_no_deform = [];

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

    if_no_deform = [if_no_deform, xlsfile{:, 22}];

end

f = figuresetting('centimeters',20,16,'times new roman',20,'on',1,'off','off');
f.figure('');
% for legends
semilogy(nan, nan, 'o','MarkerSize', 10,'MarkerEdgeColor',[229 120 148]/255,'LineWidth',2); hold on
semilogy(nan, nan, '^','MarkerSize', 9,'MarkerEdgeColor',[62 187 229]/255,'LineWidth',2); hold on
semilogy(nan, nan, 'diamond','MarkerSize', 9,'MarkerEdgeColor',[104 97 162]/255,'LineWidth',2); hold on
% plot
semilogy(FlowAngle(logical(if_C_shape)), mu_bars(logical(if_C_shape)), 'o','MarkerSize', 10,'MarkerEdgeColor',[229 120 148]/255,'LineWidth',2); hold on
semilogy(FlowAngle(logical(if_U_shape)), mu_bars(logical(if_U_shape)), 'o','MarkerSize', 10,'MarkerEdgeColor',[229 120 148]/255,'LineWidth',2); hold on
semilogy(FlowAngle(logical(if_S_shape)), mu_bars(logical(if_S_shape)), 'o','MarkerSize', 10,'MarkerEdgeColor',[229 120 148]/255,'LineWidth',2); hold on
semilogy(FlowAngle(logical(if_plus_S_shape)), mu_bars(logical(if_plus_S_shape)), 'o','MarkerSize', 10,'MarkerEdgeColor',[229 120 148]/255,'LineWidth',2); hold on
semilogy(FlowAngle(logical(if_buckled_3D)), mu_bars(logical(if_buckled_3D)), '^','MarkerSize', 9,'MarkerEdgeColor',[62 187 229]/255,'LineWidth',2); hold on
semilogy(FlowAngle(logical(if_folded)), mu_bars(logical(if_folded)), '^','MarkerSize', 9,'MarkerEdgeColor',[62 187 229]/255,'LineWidth',2); hold on
semilogy(FlowAngle(logical(if_coiled)), mu_bars(logical(if_coiled)), '^','MarkerSize', 9,'MarkerEdgeColor',[62 187 229]/255,'LineWidth',2); hold on
semilogy(FlowAngle(logical(if_no_deform)), mu_bars(logical(if_no_deform)), 'diamond','MarkerSize', 9,'MarkerEdgeColor',[104 97 162]/255,'LineWidth',2); hold on

f.interp_font('latex')
f.axes('linear',[-5 50],'log',[0 1e7],'$\alpha\,(^\circ)$','$\bar{\mu}$',24);
f.axes_ticks([0:5:45], 10.^[1:7]); 
set(gca, 'XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.2)
legend({'2D deformation','3D deformation','No deformation'}, 'FontName','Times New Roman', ...
    'Box','off','Location','southeast'); 

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Dynamics\actin_dyn_mubar.eps']);

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Dynamics\actin_dyn_mubar.pdf']);


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
semilogy(nan, nan, 'o','MarkerSize', 10,'MarkerEdgeColor',[229 120 148]/255,'LineWidth',2); hold on
semilogy(nan, nan, '^','MarkerSize', 9,'MarkerEdgeColor',[62 187 229]/255,'LineWidth',2); hold on
semilogy(nan, nan, 'diamond','MarkerSize', 9,'MarkerEdgeColor',[104 97 162]/255,'LineWidth',2); hold on
% plot
semilogy(FlowAngle(logical(if_C_shape)), mu_bars_times_flow_strength(logical(if_C_shape)), ...
    'o','MarkerSize', 10,'MarkerEdgeColor',[229 120 148]/255,'LineWidth',2); hold on
semilogy(FlowAngle(logical(if_U_shape)), mu_bars_times_flow_strength(logical(if_U_shape)), ...
    'o','MarkerSize', 10,'MarkerEdgeColor',[229 120 148]/255,'LineWidth',2); hold on
semilogy(FlowAngle(logical(if_S_shape)), mu_bars_times_flow_strength(logical(if_S_shape)), ...
    'o','MarkerSize', 10,'MarkerEdgeColor',[229 120 148]/255,'LineWidth',2); hold on
semilogy(FlowAngle(logical(if_plus_S_shape)), mu_bars_times_flow_strength(logical(if_plus_S_shape)), ...
    'o','MarkerSize', 10,'MarkerEdgeColor',[229 120 148]/255,'LineWidth',2); hold on
semilogy(FlowAngle(logical(if_buckled_3D)), mu_bars_times_flow_strength(logical(if_buckled_3D)), ...
    '^','MarkerSize', 9,'MarkerEdgeColor',[62 187 229]/255,'LineWidth',2); hold on
semilogy(FlowAngle(logical(if_folded)), mu_bars_times_flow_strength(logical(if_folded)), ...
    '^','MarkerSize', 9,'MarkerEdgeColor',[62 187 229]/255,'LineWidth',2); hold on
semilogy(FlowAngle(logical(if_coiled)), mu_bars_times_flow_strength(logical(if_coiled)), ...
    '^','MarkerSize', 9,'MarkerEdgeColor',[62 187 229]/255,'LineWidth',2); hold on
semilogy(FlowAngle(logical(if_no_deform)), mu_bars_times_flow_strength(logical(if_no_deform)), ...
    'diamond','MarkerSize', 9,'MarkerEdgeColor',[104 97 162]/255,'LineWidth',2); hold on


f.interp_font('latex')
f.axes('linear',[-5 50],'log',[100 5e7],'$\alpha\,(^\circ)$','$\bar{\mu}_{\rm m}$',24);
f.axes_ticks([0:5:45], 10.^[1:7]);
set(gca, 'XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.2)
legend({'2D deformation','3D deformation', 'No deformation'}, 'FontName','Times New Roman', ...
    'Box','off','Location','southeast'); 

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Dynamics\actin_dyn_mubar_times_delta.eps']);

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Dynamics\actin_dyn_mubar_times_delta.pdf']);


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

f = figuresetting('centimeters',20,16,'times new roman',24,'on',1,'off','off');
f.figure('');
plot(FlowAngle_passing + 2*rand(1,length(FlowAngle_passing))-1, ContourLs_passing, ...
    'o','MarkerSize', 10,'MarkerEdgeColor',[176 67 129]/255,'LineWidth',2); hold on
plot(FlowAngle_others + 2*rand(1,length(FlowAngle_others))-1, ContourLs_others, ...
    '^','MarkerSize', 9,'MarkerEdgeColor',[117 152 196]/255,'LineWidth',2); hold on
plot(FlowAngle_trapping + 2*rand(1,length(FlowAngle_trapping))-1, ContourLs_trapping, ...
    'diamond','MarkerSize', 9,'MarkerEdgeColor',[123 164 47]/255,'LineWidth', 2); hold on

f.interp_font('latex')
f.axes('linear',[-5 50],'linear',[0 80],'$\alpha\,(^\circ)$','$L\,(\rm{\mu m})$',24);
f.axes_ticks([0:5:45], [0:10:80]); grid on
set(gca,'XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.2)

legend({'Passing','Gliding','Trapping'}, 'FontName','Times New Roman', 'Box','off'); 

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Dynamics\actin_pillar_interaction_add_NOISE.pdf']);


%% Statistics about actin-pillar dynamics in PAs (mean velocity vs dyns).

clear; close all; clc;
xlsFolder = dir(['F:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Dynamics manually classification\*.xlsx']);

case_names = [];             ContourLs = [];             mu_bars = [];
FlowAngle = [];          
if_trapping = [];        if_passing = [];            if_others = [];
Ux_to_U0 = [];

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

    Ux_to_U0 = [Ux_to_U0, xlsfile{:, 19}];
 
end

% % choose an angle to plot
% chosen_angle_to_plot = 45;
% Ux_to_U0_trapping = Ux_to_U0(logical(if_trapping) & (FlowAngle == chosen_angle_to_plot));
% Ux_to_U0_passing = Ux_to_U0(logical(if_passing) & (FlowAngle == chosen_angle_to_plot));
% Ux_to_U0_others = Ux_to_U0(logical(if_others) & (FlowAngle == chosen_angle_to_plot));
% histogram(Ux_to_U0_passing, 'BinWidth', 0.3, 'Normalization', 'probability'); hold on
% histogram(Ux_to_U0_others, 'BinWidth', 0.3, 'Normalization', 'probability'); hold on
% histogram(Ux_to_U0_trapping, 'BinWidth', 0.3, 'Normalization', 'probability'); 
% legend({'Passing','Gliding','Trapping'}, 'FontName','Times New Roman', 'Box','off');

% plot all angles
Ux_to_U0_trapping = Ux_to_U0(logical(if_trapping));
Ux_to_U0_passing = Ux_to_U0(logical(if_passing));
Ux_to_U0_others = Ux_to_U0(logical(if_others));

hhh = figure('color', 'w'); set(gcf, 'Position', [100 100 800 600],'DefaultTextInterpreter', 'latex');
histogram(Ux_to_U0_passing, 'BinWidth', 0.25, 'Normalization', 'probability'); hold on
histogram(Ux_to_U0_others, 'BinWidth', 0.25, 'Normalization', 'probability'); hold on
histogram(Ux_to_U0_trapping, 'BinWidth', 0.25, 'Normalization', 'probability'); 

set(gca, 'Box', 'On', 'XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.2, ...
    'FontSize', 24, 'FontName', 'Times new roman')
legend({'Passing','Gliding','Trapping'}, 'FontName','Times New Roman', 'Box','off');

% Convert y-axis values to percentage values by multiplication
a = cellstr(num2str(get(gca,'ytick')'*100));
% Create a vector of '%' signs
pct = char(ones(size(a,1),1)*'%');
% Append the '%' signs after the percentage values
new_yticks = [char(a),pct];
for i=1:size(new_yticks,1); the_yticks{i} = new_yticks(i, :); end
% 'Reflect the changes on the plot
set(gca,'yticklabel',the_yticks)
xlabel('$\overline{U}_x / U_0$','FontSize', 24,'Interpreter', 'latex');
ylabel('Probability','FontSize', 24, 'FontName', 'Times new roman');

set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Dynamics\' ...
    'actin_pillar_interaction_vs_velocity.pdf']);



%% Statistics about different actin dynamics in PAs (partial buckling).

clear; close all; clc;
xlsFolder = dir(['F:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Dynamics manually classification\*.xlsx']);

case_names = [];             ContourLs = [];             mu_bars = [];
FlowAngle = [];          
if_partial_buckled = [];        if_global_buckled = [];   
if_has_note = [];       if_out_of_plane = [];       if_bad_reconstruction = [];
if_no_deform = [];

for ii = 1: length(xlsFolder)

    xlsfile = readcell([xlsFolder(1).folder, filesep, xlsFolder(ii).name], ...
        'Sheet','Sheet1','NumHeaderLines',1);

    case_names = [case_names; xlsfile(:, 1)];
    ContourLs = [ContourLs, xlsfile{:, 2}];
    mu_bars = [mu_bars, xlsfile{:, 14}];

    FlowAngle = [FlowAngle, xlsfile{:, 15}];

    if_has_note = [if_has_note, xlsfile{:, 11}];
    if_out_of_plane = [if_out_of_plane, xlsfile{:, 12}];
    if_bad_reconstruction = [if_bad_reconstruction, xlsfile{:, 13}];

    if_partial_buckled = [if_partial_buckled, xlsfile{:, 20}];
    if_global_buckled = [if_global_buckled, xlsfile{:, 21}];
    if_no_deform = [if_no_deform, xlsfile{:, 22}];

end

f = figuresetting('centimeters',20,16,'times new roman',20,'on',1,'off','off');
f.figure('');
% plot
plot(FlowAngle(logical(if_global_buckled)), ContourLs(logical(if_global_buckled)), ...
    'v','MarkerSize', 10,'MarkerEdgeColor',[78 101 136]/255,'LineWidth',2); hold on
plot(FlowAngle(logical(if_partial_buckled)), ContourLs(logical(if_partial_buckled)), ...
    'square','MarkerSize', 10,'MarkerEdgeColor',[227 156 109]/255,'LineWidth',2); hold on

f.interp_font('latex')
f.axes('linear',[-5 50],'linear',[0 80],'$\alpha\,(^\circ)$','$L\,(\rm{\mu m})$',24);
f.axes_ticks([0:5:45], [10:10:90]); 
set(gca, 'XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.2)
legend({'Global buckling','Partial buckling'}, 'FontName','Times New Roman', ...
    'Box','off','Location','northeast'); 

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Dynamics\actin_buckling_L.eps']);

hhh = gcf;
set(hhh,'Units','Inches');
pos = get(hhh,'Position');
set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hhh, '-dpdf',['F:\Processing & Results\' ...
    'Actin Filaments in Porous Media\Figures\Dynamics\actin_buckling_L.pdf']);



%% Statistics about actin-pillar dynamics in PAs (mean velocity vs actin length).

clear; close all; clc;
xlsFolder = dir(['F:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Dynamics manually classification\*.xlsx']);

case_names = [];             ContourLs = [];             mu_bars = [];
FlowAngle = [];          
if_trapping = [];        if_passing = [];            if_others = [];
Ux_to_U0 = [];

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

    Ux_to_U0 = [Ux_to_U0, xlsfile{:, 19}];
 
end

% % plot for different angles
for chosen_angle_to_plot = [0 10 15 20 30 35 45]
    Ux_to_U0_trapping = Ux_to_U0(logical(if_trapping) & (FlowAngle == chosen_angle_to_plot));
    Ux_to_U0_passing = Ux_to_U0(logical(if_passing) & (FlowAngle == chosen_angle_to_plot));
    Ux_to_U0_others = Ux_to_U0(logical(if_others) & (FlowAngle == chosen_angle_to_plot));

    ContourLs_trapping = ContourLs(logical(if_trapping) & (FlowAngle == chosen_angle_to_plot));
    ContourLs_passing = ContourLs(logical(if_passing) & (FlowAngle == chosen_angle_to_plot));
    ContourLs_others = ContourLs(logical(if_others) & (FlowAngle == chosen_angle_to_plot));

    hhh = figure('color', 'w'); set(gcf, 'Position', [100 100 800 600],'DefaultTextInterpreter', 'latex');
    plot(ContourLs_passing, 1./Ux_to_U0_passing, 'o','MarkerSize', 15, 'MarkerFaceColor', ...
        [176 67 129]/255, 'MarkerEdgeColor', [176 67 129]/255); hold on
    plot(ContourLs_others, 1./Ux_to_U0_others, 'o','MarkerSize', 15, 'MarkerFaceColor', ...
        [117 152 196]/255, 'MarkerEdgeColor', [117 152 196]/255); hold on
    plot(ContourLs_trapping, 1./Ux_to_U0_trapping, 'o','MarkerSize', 15, 'MarkerFaceColor', ...
        [123 164 47]/255, 'MarkerEdgeColor', [123 164 47]/255);
    hold on

    set(gca, 'Box', 'On', 'XGrid', 'On', 'YGrid', 'On', 'GridAlpha', 0.2, ...
        'FontSize', 24, 'FontName', 'Times new roman')
    legend({'Passing','Gliding','Trapping'}, 'FontName','Times New Roman', 'Box','off','Location','northwest');

    xlabel('$L\,(\rm{\mu m})$','FontSize', 24,'Interpreter', 'latex');
    ylabel('$\overline{t}_x / t_0$','FontSize', 24, 'FontName', 'Times new roman');

    set(hhh,'Units','Inches');
    pos = get(hhh,'Position');
    set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hhh, '-dpdf',['F:\Processing & Results\' ...
        'Actin Filaments in Porous Media\Figures\Dynamics\' ...
        'actin_pillar_interaction_velocity_vs_length_alpha',num2str(chosen_angle_to_plot),'deg.pdf']);
end