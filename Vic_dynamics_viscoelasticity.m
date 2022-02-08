% Classification of different morphologies which are got manually..


%% Extract the data from excel.
clear; clc; close all;
xlsfile = readcell('ForActinPostprocessing.xlsx','Sheet','Sheet2','NumHeaderLines',1); % This is the master file that contains all the information (in sheet 2).

NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.
savePath = xlsfile(:, 3);  % Saving path.

for no_Group = 1:3
    Allinfo = readtable(storePath{no_Group},'Sheet','Sheet2');  % read the 
    for ii = 5:9 % the column in the sheet
        Allinfo(:, ii) = array2table(table2array(Allinfo(:, ii))+(ii-4));
    end
    % This loop is to assign numbers to different morphologies, e.g.,
    % U-shape=1, S-shape=2, W-shape=3, Helical=4, Folded=5, Coiled=6 if they are true
    % U-shape=0, S-shape=1, W-shape=2, Helical=3, Folded=4, Coiled=5 if they are false.
    Alldata = table2array(Allinfo); % transfer to array format.
    tmp1 = find(isnan(Alldata(:, 10)));
    Alldata(tmp1', :)  = [];  % to delete the rows (or cases) which don't have the information of the contour length and so on.

    %%%%%%  different modes %%%%%%%%%%
    ToDraw_dynamicmodes = Alldata(:, [4:9, 11]);  % 4-9 are the shapes; 11 is the elastoviscous number.
    for ii = 2:6
        tmp2 = ToDraw_dynamicmodes(:, ii)<ii;
        ToDraw_dynamicmodes(tmp2', ii) = 0;
    end
    % because we assigned numbers above, here is to set the false values to
    % be 0.
    for dyn = 1:6
        tmp3 = find(ToDraw_dynamicmodes(:, dyn)==0);
        dynmodes = ToDraw_dynamicmodes(:, [dyn, 7]);  % put different modes and elastoviscous number together
        dynmodes(tmp3', :)  = [];
        % remove the false cases.
        DrawXY_dynmodes{no_Group, dyn} = dynmodes; % save the data into a cell (row: groups; column: dynamic modes).
    end

    %%%%%%  global & partial buckling %%%%%%%%%%
    ToDraw_bucklingtypes = Alldata(:, [1:2, 11]);
    for buk = 1:2
        tmp3 = find(ToDraw_bucklingtypes(:, buk)==0);
        buckletypes = ToDraw_bucklingtypes(:, [buk, 3]);  % put buckling types and elastoviscous number together
        buckletypes(tmp3', :)  = [];
        % remove the false cases.
        DrawXY_buckletypes{no_Group, buk} = buckletypes;  % save the data into a cell (row: groups; column: buckling types).
    end
end

%% Plot all points

% dynamic modes vs. elastoviscous numner
figure('color', 'w'); set(gcf, 'Position', [100 300 800 400]);
% CM1 = jet(6);  % declare the colours for each mode (6 modes)
for no_Group = 1:3
    for dyn = 1:6
        if no_Group == 3
            h1 = semilogx(DrawXY_dynmodes{no_Group, dyn}(:, 2), DrawXY_dynmodes{no_Group, dyn}(:, 1),...
                'Color', 'r', 'LineStyle', 'none', 'marker', '^', 'MarkerSize', 10); hold on
        else
            h2 = semilogx(DrawXY_dynmodes{no_Group, dyn}(:, 2), DrawXY_dynmodes{no_Group, dyn}(:, 1),...
                'Color', 'b', 'LineStyle', 'none', 'marker', 'o', 'MarkerSize', 10); hold on
        end

    end
end
ylim([0 7]);
names = {'U-Shape'; 'S-Shape'; 'W-Shape'; 'Helical'; 'Folded'; 'Coiled'};
set(gca,'ytick',[1:6],'yticklabel',names,'FontSize',10);
xlabel('$\bar{\mu}$','FontSize', 15,'Interpreter', 'latex');
ylabel('$Dynamic\ Modes$','FontSize', 14,'Interpreter', 'latex');
legend([h1 h2], '$\lambda/a=1.75$' , '$\lambda/a=1.50$', 'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,'E:\Dropbox\Research\All Plottings\General plots\dynamicmodes_elastoviscousnumber.png','Resolution',1000)

% buckling types vs. elastoviscous numner 
figure('color', 'w'); set(gcf, 'Position', [100 300 800 400]);
% CM2 = jet(2);  % declare the colours for each types (2 types)
for no_Group = 1:3
    for buk = 1:2
        if no_Group == 3
            h1 = semilogx(DrawXY_buckletypes{no_Group, buk}(:, 2), DrawXY_buckletypes{no_Group, buk}(:, 1)*buk,...
                'Color', 'r', 'LineStyle', 'none', 'marker', '^', 'MarkerSize', 10); hold on
        else
            h2 = semilogx(DrawXY_buckletypes{no_Group, buk}(:, 2), DrawXY_buckletypes{no_Group, buk}(:, 1)*buk,...
                'Color', 'b', 'LineStyle', 'none', 'marker', 'o', 'MarkerSize', 10); hold on
        end
    end
end
ylim([0 3]);
names = {'Global buckling'; 'Partial buckling'};
set(gca,'ytick',[1:2],'yticklabel',names,'FontSize',10); ytickangle(90)
xlabel('$\bar{\mu}$','FontSize', 15,'Interpreter', 'latex');
ylabel('$Dynamic\ Modes$','FontSize', 14,'Interpreter', 'latex')
legend([h1 h2], '$\lambda/a=1.75$' , '$\lambda/a=1.50$', 'Interpreter', 'latex');
% f=gcf;
% exportgraphics(f,'E:\Dropbox\Research\All Plottings\General plots\bucklingtypes_elastoviscousnumber.png','Resolution',1000)

%% Plot averaged mu_bar

% dynamic modes vs. elastoviscous numner
figure('color', 'w'); set(gcf, 'Position', [100 300 800 400]);
% CM1 = jet(6);  % declare the colours for each mode (6 modes)
for no_Group = 1:3
    for dyn = 1:6
        averaged_mu_bar_dyn(no_Group, dyn) = mean(DrawXY_dynmodes{no_Group, dyn}(:, 2));
    end
end
h1 = semilogx(averaged_mu_bar_dyn(3, :), 1:dyn,'Color', 'r', 'LineStyle', 'none', 'marker', '^', 'MarkerSize', 10); hold on
h2 = semilogx(mean(averaged_mu_bar_dyn(1:2, :), 'omitnan'), 1:dyn,'Color', 'b', 'LineStyle', 'none', 'marker', 'o', 'MarkerSize', 10); hold on
ylim([0 7]);
names = {'U-Shape'; 'S-Shape'; 'W-Shape'; 'Helical'; 'Folded'; 'Coiled'};
set(gca,'ytick',[1:6],'yticklabel',names,'FontSize',10);
xlabel('$\bar{\mu}$','FontSize', 15,'Interpreter', 'latex');
ylabel('$Dynamic\ Modes$','FontSize', 14,'Interpreter', 'latex');
legend([h1 h2], '$\lambda/a=1.75$' , '$\lambda/a=1.50$', 'Interpreter', 'latex');

% buckling types vs. elastoviscous numner 
figure('color', 'w'); set(gcf, 'Position', [100 300 800 400]);
% CM1 = jet(6);  % declare the colours for each mode (6 modes)
for no_Group = 1:3
    for buk = 1:2
        averaged_mu_bar_buk(no_Group, buk) = mean(DrawXY_buckletypes{no_Group, buk}(:, 2));
    end
end
h1 = semilogx(averaged_mu_bar_buk(3, :), 1:buk,'Color', 'r', 'LineStyle', 'none', 'marker', '^', 'MarkerSize', 10); hold on
h2 = semilogx(mean(averaged_mu_bar_buk(1:2, :), 'omitnan'), 1:buk,'Color', 'b', 'LineStyle', 'none', 'marker', 'o', 'MarkerSize', 10); hold on
ylim([0 7]);
names = {'U-Shape'; 'S-Shape'; 'W-Shape'; 'Helical'; 'Folded'; 'Coiled'};
set(gca,'ytick',[1:6],'yticklabel',names,'FontSize',10);
xlabel('$\bar{\mu}$','FontSize', 15,'Interpreter', 'latex');
ylabel('$Dynamic\ Modes$','FontSize', 14,'Interpreter', 'latex');
legend([h1 h2], '$\lambda/a=1.75$' , '$\lambda/a=1.50$', 'Interpreter', 'latex');