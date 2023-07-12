%% Plot the bending energy, Lee, Chi, and Aniso together with normalized CoM_x positions.

clear; close all; clc;

set(0,'DefaultAxesFontSize',14);
set(0,'defaulttextfontsize',14);
set(0,'defaultAxesFontName', 'times new roman');
set(0,'defaultTextFontName', 'times new roman');
set(0,'defaulttextInterpreter','latex') 

mag = 0.1; % um/pixel
cmap = cmocean('amp'); 

xlsfile = readcell('ForActinPostprocessing.xlsx','Sheet','Sheet1','NumHeaderLines',1);
% This is the file that contains all the information about the later processing (in sheet 1).

NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.
Array_angles = cell2mat(xlsfile(:, 14));  % The flow angles.

for current_angle = [0 10 15 20 30 35 45] % different angles

% %     Plot_CoM_x_all = []; Plot_Chi_all = [];
    x_min = 0; x_max = 0;

    fiber_CoM_xy_in_lattice_all = [];
    figure('color', 'w'); set(gcf, 'Position', [100 100 900 400]);

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

        excelpathname = 'F:\Processing & Results\Actin Filaments in Porous Media\Dynamics manually classification\';
        excelname = ['Dynamics ',num2str(the_exp_date),'-Actin.xlsx']; % read excel
        xlsfile = readcell(fullfile(excelpathname, excelname),'Sheet','Sheet1','NumHeaderLines',1);
        filenames_xls = xlsfile(:, 1);

        for file_ind = 1:length(thefiles) % load the cases

            filename = thefiles(file_ind).name;

            if contains(filename, 'PlusInfo_')

                load(fullfile(thefiles(1).folder, thefiles(file_ind).name));

                name_ind = find(cellfun(@(x) contains(thefiles(file_ind).name, x), filenames_xls)); % index in the excel

                if xlsfile{name_ind, 12} || xlsfile{name_ind, 13}

                    continue

                else

                    spl_Ls = xy.arclen_spl(Good_case_frm);
                    L_0 = VicFc_Get_ContourLength(spl_Ls) * mag; % get the filament length

                    if round(L_0 * 5 - 10) > 256
                        color_ind = 256;
                    else
                        color_ind = round(L_0 * 5 - 10);
                    end

                    P_inter = InterX([0 2000; 0 0], [Rotated_fiber_CoM_xy_normalized(1, :); ...
                        Rotated_fiber_CoM_xy_normalized(2, :)]);

                    if ~isempty(P_inter)

                        % plot based on the pillar array coordinates (rotated back to unit cell)
                        Plot_CoM_x = Rotated_fiber_CoM_xy_normalized(1, :)' - P_inter(1);
                        plot((Plot_CoM_x+lzero)*mag, Chi/pi,'.', 'LineStyle','none', 'MarkerSize', 10, 'Color', cmap(color_ind, :));
                        hold on
                        x_min = min( min((Plot_CoM_x+lzero)*mag), x_min );
                        x_max = max( max((Plot_CoM_x+lzero)*mag), x_max );

%                         Rotated_fiber_CoM_xy_normalized(1, :) = Rotated_fiber_CoM_xy_normalized(1, :) - P_inter(1);
%                         fiber_CoM_xy_normalized = RotMatrix_correct \ Rotated_fiber_CoM_xy_normalized;
%                         plot((fiber_CoM_xy_normalized(1, :)'+lzero)*mag, Chi/pi,'.', 'LineStyle','none', 'MarkerSize', 10, 'Color', cmap(color_ind, :));
%                         hold on

% %                         Plot_CoM_x_all = [(Plot_CoM_x+lzero)*mag; Plot_CoM_x_all];
% %                         Plot_Chi_all = [Chi/pi; Plot_Chi_all];
                    else
                        Plot_CoM_x = Rotated_fiber_CoM_xy_normalized(1, :)';
                        plot((Plot_CoM_x+lzero)*mag, Chi/pi,'.', 'LineStyle','none', 'MarkerSize', 10, 'Color', cmap(color_ind, :));
                        hold on
                        x_min = min( min((Plot_CoM_x+lzero)*mag), x_min );
                        x_max = max( max((Plot_CoM_x+lzero)*mag), x_max );

%                         fiber_CoM_xy_normalized = RotMatrix_correct \ Rotated_fiber_CoM_xy_normalized;
%                         plot((fiber_CoM_xy_normalized(1, :)'+lzero)*mag, Chi/pi,'.', 'LineStyle','none', 'MarkerSize', 10, 'Color', cmap(color_ind, :));
%                         hold on

% %                         Plot_CoM_x_all = [(Plot_CoM_x+lzero)*mag; Plot_CoM_x_all];
% %                         Plot_Chi_all = [Chi/pi; Plot_Chi_all];
                    end

                    clearvars P_inter
                end

% % %                 D = pdist2([0 0; 0 1; 1 0; 1 1], fiber_CoM_xy_in_lattice);
% % % 
% % %                 if min(min(D)) > mean(radii)/mean([Ctr2Ctr_x, Ctr2Ctr_y])
% % %                     Plot_CoM_x = Rotated_fiber_CoM_xy_normalized(1, :)';
% % %                     plot((Plot_CoM_x+lzero)*mag, Chi/pi,'.', 'LineStyle','none', 'MarkerSize', 10);
% % %                     hold on
% % %                 end


            end
        end
    end
    ylim([-0.5 0.5]); xlim([x_min-10 x_max+10]);
    ylabel('$\theta / \pi$', 'FontSize', 24);
    xlabel('$x_{\rm array}\,(\rm{\mu m})$', 'FontSize', 24);
    for jj = -7:8
        plot([Ctr2Ctr_x*jj*mag Ctr2Ctr_x*jj*mag], [-0.5 0.5], ':k','LineWidth',1); hold on
    end
    text(x_min-7, 0.39, ['$\alpha =',num2str(current_angle),'^{\circ}$'],'FontSize', 24, ...
    'Interpreter', 'latex','BackgroundColor',[.7 .7 .7])
    set(gca,'Box', 'On','FontSize', 24)

    if current_angle == 45
        caxis([0 50]); cmocean('amp');
        hcb=colorbar;
        hcb.Label.String = '$L\ (\mathrm{\mu m})$';
        hcb.Label.Interpreter = 'LaTeX';
        hcb.TickLabelInterpreter = 'LaTeX';
        hcb.FontSize = 24;
        hcb.Location = 'eastoutside';
    end

    hhh = gcf;
    set(hhh,'Units','Inches');
    pos = get(hhh,'Position');
    set(hhh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(hhh, '-dpdf',['F:\Processing & Results\Actin Filaments in Porous Media' ...
        '\Figures\Periodic plots\Fiber_ori_theta_TiltAngle=',num2str(current_angle),'.pdf']);

%%%%%%%%% bin and plot all data together
% %     % sort the results in ascending order
% %     [Plot_CoM_x_all,index_CoM_x] = sort(Plot_CoM_x_all);
% %     Plot_Chi_all = Plot_Chi_all(index_CoM_x);
% % 
% %     % data binning
% %     pitch = 2; % mean pitch of the bin
% %     edges = min(Plot_CoM_x_all)-1/2*pitch : pitch : max(Plot_CoM_x_all)+1/2*pitch; % array of bin where to count
% %     subs_array = discretize(Plot_CoM_x_all,edges); % discretize the data such as every bin in Plot_CoM_x_all has an identical number
% % 
% %     % results of the binning operation
% %     try
% %         Plot_CoM_x_all_in_bin = splitapply(@mean,Plot_CoM_x_all,subs_array);
% %         Plot_Chi_all_in_bin = splitapply(@mean,Plot_Chi_all,subs_array);
% %         Plot_CoM_x_all_in_bin_std = splitapply(@std,Plot_CoM_x_all,subs_array);
% %         Plot_Chi_all_in_bin_std = splitapply(@std,Plot_Chi_all,subs_array);
% % 
% %     catch
% %         subs_array = subs_array(1: round(0.97* numel(subs_array)));
% %         Plot_Chi_all = Plot_Chi_all(1: round(0.97* numel(Plot_Chi_all)));
% %         Plot_CoM_x_all = Plot_CoM_x_all(1: round(0.97* numel(Plot_CoM_x_all)));
% % 
% %         Plot_CoM_x_all_in_bin = splitapply(@mean,Plot_CoM_x_all,subs_array);
% %         Plot_Chi_all_in_bin = splitapply(@mean,Plot_Chi_all,subs_array);
% %         Plot_CoM_x_all_in_bin_std = splitapply(@std,Plot_CoM_x_all,subs_array);
% %         Plot_Chi_all_in_bin_std = splitapply(@std,Plot_Chi_all,subs_array);
% %     end
% % 
% %     figure; plot(Plot_CoM_x_all_in_bin, Plot_Chi_all_in_bin)


end
