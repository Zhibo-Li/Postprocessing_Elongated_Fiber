clear; close all; clc

%%%
scalefactor = 0.1;  % um/pix
delta_t = 0.02; % second

% load the data
data_path = uigetdir('F:\Processing & Results\FSI - Actin &  Individual Obstacle\', ...
    'Choose the data path which contains the *.mat results !');
Files = dir([data_path,'\*.mat']);

theimgs = dir(['F:\Experimental Data (EXTRACTED)\FSI - Actin &  Individual Obstacle' ...
    '\',Files(1).folder(60:end-8),'\AfterAveBGR\*.tif']);

for file_ind = 1:length(Files)

    filename = Files(file_ind).name;

    if contains(filename, 'AddInfo_')

        load(fullfile(Files(1).folder, Files(file_ind).name));
        obs_x = mean(obs_2d(:,2)) * scalefactor;

        Good_case_frm = find(ismember(xy(1).frame, Good_case));
        for ii = 1: size(Good_case,2)-1

            velocity_ins(ii) = (xy(1).centroid{Good_case_frm(ii+1)}(1) - xy(1).centroid{Good_case_frm(ii)}(1)) ...
                / (Good_case(ii+1) - Good_case(ii)) * scalefactor / delta_t;
            pos_x(ii) = (xy(1).centroid{Good_case_frm(ii+1)}(1) + xy(1).centroid{Good_case_frm(ii)}(1)) / 2 * scalefactor;

        end
        figure; plot(pos_x, movmean(velocity_ins, 5), 'ro', 'MarkerSize', 7);
        xlabel('$X\ (\mu{m})$','Interpreter', 'latex');
        ylabel('$V_x\ (\mu{m}/s)$','Interpreter', 'latex');
        xlim([0 200]); ylim([0 300])

        hold on
        line([obs_x obs_x], [0 1000], 'color', 'k', 'linewidth', 3)
        annotation('textarrow', [0.64 0.54], [0.81 0.86],'String','$X_{\rm obs,mid}$', ...
            'FontSize', 24, 'Interpreter', 'latex')

        set_plot(gcf, gca)

        f=gcf;
        exportgraphics(f,['F:\Processing & Results\FSI - Actin &  Individual Obstacle\' ...
            '20230330-Actin-Individual_triangularPillar_uppoint\velocity_X-position\' ...
            , Files(case_No).name(1: end-17),'.png'],'Resolution',100)

        close
        clearvars pos_x velocity_ins

    end
end