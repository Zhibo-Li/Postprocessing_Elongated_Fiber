% Classification of different morphologies which are got manually..
% No old mask data is included. (20220417)

%% Extract the data from excel.
clear; clc; close all;
xlsfile = readcell('ForActinPostprocessing.xlsx','Sheet','Sheet2','NumHeaderLines',1); 
%%% This is the master file that contains all the information (in sheet 2).

NumGroup = size(xlsfile, 1)-2;  % Number of the groups to be calculated.
ExpDate = xlsfile(3:end, 1);  % The experiment date.
storePath = xlsfile(3:end, 2);  % Path of the data to be processed.
savePath = xlsfile(3:end, 3);  % Saving path.
PAtype = xlsfile(3:end, 4); % Pillar array types.
PAtype_label = xlsfile(3:end, 5); % PAtype labels (in number). Easier for drawing.
PAtype_label = cell2mat(PAtype_label)-2; % convert format.

CM = jet(max(PAtype_label));  % declare the colours for different pillar array types.
MyMarkers = {'p','o','x','s','d','v','>','<','+','h','*','.','^'};  % declare the markers for different pillar array types.

% Plot histogram of contour lengths in different PA types.
ContourL_all = cell(1, max(PAtype_label));
PAtype_name = cell(1, max(PAtype_label));
for no_Group = 1: NumGroup
    Allinfo = readtable(storePath{no_Group},'Sheet','Sheet2','VariableNamingRule','preserve');  % read the 
    ContourL = table2array(Allinfo(:, 10)); % transfer to array format (contour length).  
    ContourL(isnan(ContourL)) = []; % remove the NAN results.
    ContourL(ContourL<4) = []; ContourL(ContourL>50) = []; % Select cases that contour length between 10um and 24um.
    type_ind = PAtype_label(no_Group); % PA type indicators.
    ContourL_all(type_ind) = {[ContourL_all{type_ind};ContourL]}; % classify countour according to different PA types.
    PAtype_name(type_ind) = PAtype(no_Group);
end

% Plot individual PDFs.
% for PAtype_ind = 1: max(PAtype_label)
%     figure('color', 'w'); set(gcf, 'Position', [100 300 800 400]);
%     histogram(ContourL_all{PAtype_ind}, 'BinWidth',2,'BinLimits',[4,50]); 
%     title(PAtype_name{PAtype_ind}, Interpreter='latex',FontSize=12);
%     set(gca,'FontSize',12);
%     xlabel('$Contour\ length\ ({\mu}m)$','FontSize', 12,'Interpreter', 'latex');
%     ylabel('$Number\ of\ events$','FontSize', 12,'Interpreter', 'latex');
% %     f=gcf;
% %     exportgraphics(f,['E:\Dropbox\Research\All Plottings\General plots\',PAtype_name{PAtype_ind},'_contourL_histogram.png'],'Resolution',500)
% end

% Plot PDF together.
edges = 4:4:52;
for PAtype_ind = 1: max(PAtype_label)
    data = ContourL_all{PAtype_ind};
    h{PAtype_ind} = histcounts(data,edges);
end 
figure('color', 'w'); set(gcf, 'Position', get(0, 'Screensize'));
bar(edges(1:end-1),[h{1}; h{2}; h{3}; h{4};h{5}]');
set(gca,'FontSize',18);
xlabel('$Contour\ length\ ({\mu}m)$','FontSize', 22,'Interpreter', 'latex');
ylabel('$Number\ of\ events$','FontSize', 22,'Interpreter', 'latex');
legend('$Rhombic\ \lambda/a=1.50$', '$Rhombic\ \lambda/a=1.75$','$Square\ \lambda/a=1.50\ Flow\ angle=0^{\circ}$',...
    '$Square\ \lambda/a=1.50\ Flow\ angle=10^{\circ}$','$Square\ \lambda/a=1.50\ Flow\ angle=20^{\circ}$','FontSize', 18,'Interpreter', 'latex', 'Location','best');  
% f=gcf;
% exportgraphics(f,'D:\Dropbox\Research\All Plottings\General plots\pdf_of_cases.png','Resolution',100)



%% Extract data for the following plotting.
for no_Group = 1: NumGroup
    Allinfo = readtable(storePath{no_Group},'Sheet','Sheet2','VariableNamingRule','preserve');  % read the 
    for ii = 5:9 % the column in the sheet
        Allinfo(:, ii) = array2table(table2array(Allinfo(:, ii))+(ii-4));
    end
    % This loop is to assign numbers to different morphologies, e.g.,
    % U-shape=1, S-shape=2, W-shape=3, Helical=4, Folded=5, Coiled=6 if they are true
    % U-shape=0, S-shape=1, W-shape=2, Helical=3, Folded=4, Coiled=5 if they are false.
    Alldata = table2array(Allinfo(:, 1:11)); % transfer to array format (keep only useful data).
    tmp1 = find(isnan(Alldata(:, 10)));  % Note: all the 'tmp's in this code is used as matrix index.
    Alldata(tmp1', :)  = [];  % to delete the rows (or cases) which don't have the information of the contour length and so on.

    %%%%%%  different modes %%%%%%%%%%
    ToDraw_dynamicmodes = Alldata(:, [4:9, 10:11]);  % 4-9 are the shapes; 11 is the elastoviscous number.
    ToDraw_dynamicmodes(ToDraw_dynamicmodes(:,7)<14, :) = []; % Select cases that contour length between 10um and 24um. 7th column is contour length.
    ToDraw_dynamicmodes(ToDraw_dynamicmodes(:,7)>38, :) = []; % Select cases that contour length between 10um and 24um.
    for ii = 2:6
        tmp2 = ToDraw_dynamicmodes(:, ii)<ii;
        ToDraw_dynamicmodes(tmp2', ii) = 0;
    end
    % because we assigned numbers above, here is to set the false values to
    % be 0.
    for dyn = 1:6
        tmp3 = find(ToDraw_dynamicmodes(:, dyn)==0);
        dynmodes = ToDraw_dynamicmodes(:, [dyn, 8]);  % put different modes and contour length (7) or elastoviscous number (8) together.
        dynmodes(tmp3', :)  = [];
        % remove the false cases.
        DrawXY_dynmodes{no_Group, dyn} = dynmodes; % save the data into a cell (row: groups; column: dynamic modes).
    end
    DrawXY_dynmodes{no_Group, dyn+1} = PAtype_label(no_Group);

    %%%%%%  global & partial buckling %%%%%%%%%%
    ToDraw_bucklingtypes = Alldata(:, [1:2, 10:11]);
    for buk = 1:2
        tmp3 = find(ToDraw_bucklingtypes(:, buk)==0);
        buckletypes = ToDraw_bucklingtypes(:, [buk, 4]);  % put buckling types and contour length (3) or elastoviscous number (4) together.
        buckletypes(tmp3', :)  = [];
        % remove the false cases.
        DrawXY_buckletypes{no_Group, buk} = buckletypes;  % save the data into a cell (row: groups; column: buckling types).
    end
    DrawXY_buckletypes{no_Group, buk+1} = PAtype_label(no_Group);
end


%% Plot all points
% dynamic modes vs. elastoviscous numner
figure('color', 'w'); set(gcf, 'Position', [100 300 800 400]);
for no_Group = 1:NumGroup
    PAtype_labels_ind = DrawXY_dynmodes{no_Group, 7};
    for dyn = 1:6
        if ~isempty(DrawXY_dynmodes{no_Group, dyn}) 
            h_dyn(PAtype_labels_ind) = semilogx(DrawXY_dynmodes{no_Group, dyn}(:, 2), ...
                DrawXY_dynmodes{no_Group, dyn}(:, 1), 'Color', CM(PAtype_labels_ind,:), ...
                'LineStyle', 'none', 'marker', MyMarkers{1, PAtype_labels_ind}, 'MarkerSize', 10,'LineWidth', 2); hold on
        end
    end
end
ylim([0 7]);
names = {'U-Shape'; 'S-Shape'; 'W-Shape'; 'Helical'; 'Folded'; 'Coiled'};
set(gca,'ytick',[1:6],'yticklabel',names,'FontSize',10);
xlabel('$\bar{\mu}$','FontSize', 15,'Interpreter', 'latex');
ylabel('$Dynamic\ Modes$','FontSize', 14,'Interpreter', 'latex');
legend(h_dyn,'$Rhombic\ \lambda/a=1.50$', '$Rhombic\ \lambda/a=1.75$','$Square\ \lambda/a=1.50\ Flow\ angle=0^{\circ}$',...
    '$Square\ \lambda/a=1.50\ Flow\ angle=10^{\circ}$','$Square\ \lambda/a=1.50\ Flow\ angle=20^{\circ}$','Interpreter', 'latex', 'Location','best');  
%%% change the legends according to the excel 'ForActinPostprocessing.xlsx - sheet2'.
% f=gcf;
% exportgraphics(f,'E:\Dropbox\Research\All Plottings\General plots\dynamicmodes_elastoviscousnumber.png','Resolution',1000)

% buckling types vs. elastoviscous numner 
figure('color', 'w'); set(gcf, 'Position', [100 300 800 400]);
for no_Group = 1:NumGroup
    PAtype_labels_ind = DrawXY_buckletypes{no_Group, 3};
    for buk = 1:2
        if ~isempty(DrawXY_buckletypes{no_Group, buk}) 
            h_buk(PAtype_labels_ind) = semilogx(DrawXY_buckletypes{no_Group, buk}(:, 2), ...
                DrawXY_buckletypes{no_Group, buk}(:, 1)*buk, 'Color', CM(PAtype_labels_ind,:), ...
                'LineStyle', 'none', 'marker', MyMarkers{1, PAtype_labels_ind}, 'MarkerSize', 10,'LineWidth', 2);
            hold on  % *buk: global buckling (buk=1); partial buckling (buk=2).
        end
    end
end
ylim([0 3]);
names = {'Global buckling'; 'Partial buckling'};
set(gca,'ytick',[1:2],'yticklabel',names,'FontSize',10); ytickangle(90)
xlabel('$\bar{\mu}$','FontSize', 15,'Interpreter', 'latex');
ylabel('$Dynamic\ Modes$','FontSize', 14,'Interpreter', 'latex');
legend(h_buk, '$Rhombic\ \lambda/a=1.50$', '$Rhombic\ \lambda/a=1.75$','$Square\ \lambda/a=1.50\ Flow\ angle=0^{\circ}$',...
    '$Square\ \lambda/a=1.50\ Flow\ angle=10^{\circ}$','$Square\ \lambda/a=1.50\ Flow\ angle=20^{\circ}$','Interpreter','latex', 'Location','best');  
%%% change the legends according to the excel 'ForActinPostprocessing.xlsx - sheet2'.
% f=gcf;
% exportgraphics(f,'E:\Dropbox\Research\All Plottings\General plots\bucklingtypes_elastoviscousnumber.png','Resolution',1000)

%% Plot averaged mu_bar

% dynamic modes vs. elastoviscous numner
figure('color', 'w'); set(gcf, 'Position', [100 300 800 400]);
for no_Group = 1:NumGroup
    for dyn = 1:6
        averaged_mu_bar_dyn(no_Group, dyn) = mean(DrawXY_dynmodes{no_Group, dyn}(:, 2), 'omitnan');
    end
end
for PAtypeloop = 1:max(PAtype_label)  % Loop to draw according to the PA type labels.
    tmp4 = find(PAtype_label==PAtypeloop);
    Todraw_averaged_mu_bar_dyn(PAtypeloop, :) = mean(averaged_mu_bar_dyn(tmp4, :),1);
    h_dyn_avg(PAtypeloop) = semilogx(Todraw_averaged_mu_bar_dyn(PAtypeloop, :), 1:dyn,'Color', CM(PAtypeloop,:), ...
        'LineStyle', 'none', 'marker', MyMarkers{1, PAtypeloop}, 'MarkerSize', 10,'LineWidth', 2); hold on
end
ylim([0 7]);
names = {'U-Shape'; 'S-Shape'; 'W-Shape'; 'Helical'; 'Folded'; 'Coiled'};
set(gca,'ytick',[1:6],'yticklabel',names,'FontSize',10);
xlabel('$\bar{\mu}$','FontSize', 15,'Interpreter', 'latex');
ylabel('$Dynamic\ Modes$','FontSize', 14,'Interpreter', 'latex');
legend(h_dyn_avg, '$Rhombic\ \lambda/a=1.50$', '$Rhombic\ \lambda/a=1.75$','$Square\ \lambda/a=1.50\ Flow\ angle=0^{\circ}$',...
    '$Square\ \lambda/a=1.50\ Flow\ angle=10^{\circ}$','$Square\ \lambda/a=1.50\ Flow\ angle=20^{\circ}$','Interpreter', 'latex', 'Location','best');  
%%% change the legends according to the excel 'ForActinPostprocessing.xlsx - sheet2'.
% f=gcf;
% exportgraphics(f,'E:\Dropbox\Research\All Plottings\General plots\dynamicmodes_AVG-elastoviscousnumber_add20220216data.png','Resolution',500)

% buckling types vs. elastoviscous numner 
figure('color', 'w'); set(gcf, 'Position', [100 300 800 400]);
for no_Group = 1:NumGroup
    for buk = 1:2
        averaged_mu_bar_buk(no_Group, buk) = mean(DrawXY_buckletypes{no_Group, buk}(:, 2), 'omitnan');
    end
end
for PAtypeloop = 1:max(PAtype_label)  % Loop to draw according to the PA type labels.
    tmp5 = find(PAtype_label==PAtypeloop);
    Todraw_averaged_mu_bar_buk(PAtypeloop, :) = mean(averaged_mu_bar_buk(tmp5, :),1);
    h_buk_avg(PAtypeloop) = semilogx(Todraw_averaged_mu_bar_buk(PAtypeloop, :), 1:buk,'Color', CM(PAtypeloop,:), ...
        'LineStyle', 'none', 'marker', MyMarkers{1, PAtypeloop}, 'MarkerSize', 10,'LineWidth', 2); hold on
end
ylim([0 3]);
names = {'Global buckling'; 'Partial buckling'};
set(gca,'ytick',[1:2],'yticklabel',names,'FontSize',10);ytickangle(90)
xlabel('$\bar{\mu}$','FontSize', 15,'Interpreter', 'latex');
ylabel('$Dynamic\ Modes$','FontSize', 14,'Interpreter', 'latex');
legend(h_buk_avg, '$Rhombic\ \lambda/a=1.50$', '$Rhombic\ \lambda/a=1.75$','$Square\ \lambda/a=1.50\ Flow\ angle=0^{\circ}$',...
    '$Square\ \lambda/a=1.50\ Flow\ angle=10^{\circ}$','$Square\ \lambda/a=1.50\ Flow\ angle=20^{\circ}$','Interpreter', 'latex', 'Location','best');  
%%% change the legends according to the excel 'ForActinPostprocessing.xlsx - sheet2'.
% f=gcf;
% exportgraphics(f,'E:\Dropbox\Research\All Plottings\General plots\bucklingtypes_AVG-elastoviscousnumber_add20220216data.png','Resolution',500)













%% Old version
% Extract the data from excel.
clear; clc; close all;
xlsfile = readcell('ForActinPostprocessing.xlsx','Sheet','Sheet2','NumHeaderLines',1); 
%%% This is the master file that contains all the information (in sheet 2).

NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.
savePath = xlsfile(:, 3);  % Saving path.
PAtype = xlsfile(:, 4); % Pillar array types.
PAtype_label = xlsfile(:, 5); % PAtype labels (in number). Easier for drawing.
PAtype_label = cell2mat(PAtype_label); % convert format.


%% Plot histogram of contour lengths in different PA types.
ContourL_all = cell(1, max(PAtype_label));
PAtype_name = cell(1, max(PAtype_label));
for no_Group = 1: NumGroup
    Allinfo = readtable(storePath{no_Group},'Sheet','Sheet2','VariableNamingRule','preserve');  % read the 
    ContourL = table2array(Allinfo(:, 10)); % transfer to array format (contour length).  
    ContourL(isnan(ContourL))  = []; % remove the NAN results.
    type_ind = PAtype_label(no_Group); % PA type indicators.
    ContourL_all(type_ind) = {[ContourL_all{type_ind};ContourL]}; % classify countour according to different PA types.
    PAtype_name(type_ind) = PAtype(no_Group);
end

for PAtype_ind = 1: max(PAtype_label)
    figure('color', 'w'); set(gcf, 'Position', [100 300 800 400]);
    histogram(ContourL_all{PAtype_ind}, 5); 
    title(PAtype_name{PAtype_ind}, Interpreter='latex',FontSize=12);
    set(gca,'FontSize',12);
    xlabel('$Contour\ length\ ({\mu}m)$','FontSize', 12,'Interpreter', 'latex');
    ylabel('$Number\ of\ events$','FontSize', 12,'Interpreter', 'latex');
%     f=gcf;
%     exportgraphics(f,['E:\Dropbox\Research\All Plottings\General plots\',PAtype_name{PAtype_ind},'_contourL_histogram.png'],'Resolution',500)
end


%% Extract data for the following plotting.
for no_Group = 1: NumGroup
    Allinfo = readtable(storePath{no_Group},'Sheet','Sheet2','VariableNamingRule','preserve');  % read the 
    for ii = 5:9 % the column in the sheet
        Allinfo(:, ii) = array2table(table2array(Allinfo(:, ii))+(ii-4));
    end
    % This loop is to assign numbers to different morphologies, e.g.,
    % U-shape=1, S-shape=2, W-shape=3, Helical=4, Folded=5, Coiled=6 if they are true
    % U-shape=0, S-shape=1, W-shape=2, Helical=3, Folded=4, Coiled=5 if they are false.
    Alldata = table2array(Allinfo(:, 1:11)); % transfer to array format (keep only useful data).
    tmp1 = find(isnan(Alldata(:, 10)));  % Note: all the 'tmp's in this code is used as matrix index.
    Alldata(tmp1', :)  = [];  % to delete the rows (or cases) which don't have the information of the contour length and so on.

    %%%%%%  different modes %%%%%%%%%%
    ToDraw_dynamicmodes = Alldata(:, [4:9, 10]);  % 4-9 are the shapes; 11 is the elastoviscous number.
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
    DrawXY_dynmodes{no_Group, dyn+1} = PAtype_label(no_Group);

    %%%%%%  global & partial buckling %%%%%%%%%%
    ToDraw_bucklingtypes = Alldata(:, [1:2, 10]);
    for buk = 1:2
        tmp3 = find(ToDraw_bucklingtypes(:, buk)==0);
        buckletypes = ToDraw_bucklingtypes(:, [buk, 3]);  % put buckling types and elastoviscous number together
        buckletypes(tmp3', :)  = [];
        % remove the false cases.
        DrawXY_buckletypes{no_Group, buk} = buckletypes;  % save the data into a cell (row: groups; column: buckling types).
    end
    DrawXY_buckletypes{no_Group, buk+1} = PAtype_label(no_Group);
end


%% Plot all points
% dynamic modes vs. elastoviscous numner
figure('color', 'w'); set(gcf, 'Position', [100 300 800 400]);
CM = jet(max(PAtype_label));  % declare the colours for different pillar array types.
MyMarkers = {'p','o','x','s','d','v','>','<','+','h','*','.','^'};  % declare the markers for different pillar array types.
for no_Group = 1:NumGroup
    PAtype_labels_ind = DrawXY_dynmodes{no_Group, 7};
    for dyn = 1:6
        if ~isempty(DrawXY_dynmodes{no_Group, dyn}) 
            h_dyn(PAtype_labels_ind) = semilogx(DrawXY_dynmodes{no_Group, dyn}(:, 2), ...
                DrawXY_dynmodes{no_Group, dyn}(:, 1), 'Color', CM(PAtype_labels_ind,:), ...
                'LineStyle', 'none', 'marker', MyMarkers{1, PAtype_labels_ind}, 'MarkerSize', 10,'LineWidth', 2); hold on
        end
    end
end
ylim([0 7]);
names = {'U-Shape'; 'S-Shape'; 'W-Shape'; 'Helical'; 'Folded'; 'Coiled'};
set(gca,'ytick',[1:6],'yticklabel',names,'FontSize',10);
xlabel('$\bar{\mu}$','FontSize', 15,'Interpreter', 'latex');
ylabel('$Dynamic\ Modes$','FontSize', 14,'Interpreter', 'latex');
legend(h_dyn, '$Square\ \lambda/a=1.88$' , '$Square\ \lambda/a=1.56$',...
    '$Rhombic\ \lambda/a=1.50$', '$Rhombic\ \lambda/a=1.75$','$Square\ \lambda/a=1.50$',...
    '$Square\ \lambda/a=1.50\ Flow\ angle=10^{\circ}$','Interpreter', 'latex', 'Location','best');  
%%% change the legends according to the excel 'ForActinPostprocessing.xlsx - sheet2'.
% f=gcf;
% exportgraphics(f,'E:\Dropbox\Research\All Plottings\General plots\dynamicmodes_elastoviscousnumber.png','Resolution',1000)

% buckling types vs. elastoviscous numner 
figure('color', 'w'); set(gcf, 'Position', [100 300 800 400]);
for no_Group = 1:NumGroup
    PAtype_labels_ind = DrawXY_buckletypes{no_Group, 3};
    for buk = 1:2
        if ~isempty(DrawXY_buckletypes{no_Group, buk}) 
            h_buk(PAtype_labels_ind) = semilogx(DrawXY_buckletypes{no_Group, buk}(:, 2), ...
                DrawXY_buckletypes{no_Group, buk}(:, 1)*buk, 'Color', CM(PAtype_labels_ind,:), ...
                'LineStyle', 'none', 'marker', MyMarkers{1, PAtype_labels_ind}, 'MarkerSize', 10,'LineWidth', 2);
            hold on  % *buk: global buckling (buk=1); partial buckling (buk=2).
        end
    end
end
ylim([0 3]);
names = {'Global buckling'; 'Partial buckling'};
set(gca,'ytick',[1:2],'yticklabel',names,'FontSize',10); ytickangle(90)
xlabel('$\bar{\mu}$','FontSize', 15,'Interpreter', 'latex');
ylabel('$Dynamic\ Modes$','FontSize', 14,'Interpreter', 'latex');
legend(h_buk, '$Square\ \lambda/a=1.88$' , '$Square\ \lambda/a=1.56$',...
    '$Rhombic\ \lambda/a=1.50$', '$Rhombic\ \lambda/a=1.75$','$Square\ \lambda/a=1.50$',...
    '$Square\ \lambda/a=1.50\ Flow\ angle=10^{\circ}$','Interpreter','latex', 'Location','best');  
%%% change the legends according to the excel 'ForActinPostprocessing.xlsx - sheet2'.
% f=gcf;
% exportgraphics(f,'E:\Dropbox\Research\All Plottings\General plots\bucklingtypes_elastoviscousnumber.png','Resolution',1000)

%% Plot averaged mu_bar

% dynamic modes vs. elastoviscous numner
figure('color', 'w'); set(gcf, 'Position', [100 300 800 400]);
for no_Group = 1:NumGroup
    for dyn = 1:6
        averaged_mu_bar_dyn(no_Group, dyn) = mean(DrawXY_dynmodes{no_Group, dyn}(:, 2), 'omitnan');
    end
end
for PAtypeloop = 1:max(PAtype_label)  % Loop to draw according to the PA type labels.
    tmp4 = find(PAtype_label==PAtypeloop);
    Todraw_averaged_mu_bar_dyn(PAtypeloop, :) = mean(averaged_mu_bar_dyn(tmp4, :),1);
    h_dyn_avg(PAtypeloop) = semilogx(Todraw_averaged_mu_bar_dyn(PAtypeloop, :), 1:dyn,'Color', CM(PAtypeloop,:), ...
        'LineStyle', 'none', 'marker', MyMarkers{1, PAtypeloop}, 'MarkerSize', 10,'LineWidth', 2); hold on
end
ylim([0 7]);
names = {'U-Shape'; 'S-Shape'; 'W-Shape'; 'Helical'; 'Folded'; 'Coiled'};
set(gca,'ytick',[1:6],'yticklabel',names,'FontSize',10);
xlabel('$\bar{\mu}$','FontSize', 15,'Interpreter', 'latex');
ylabel('$Dynamic\ Modes$','FontSize', 14,'Interpreter', 'latex');
legend(h_dyn_avg, '$Square\ \lambda/a=1.88$' , '$Square\ \lambda/a=1.56$',...
    '$Rhombic\ \lambda/a=1.50$', '$Rhombic\ \lambda/a=1.75$','$Square\ \lambda/a=1.50$',...
    '$Square\ \lambda/a=1.50\ Flow\ angle=10^{\circ}$','Interpreter', 'latex', 'Location','best');  
%%% change the legends according to the excel 'ForActinPostprocessing.xlsx - sheet2'.
% f=gcf;
% exportgraphics(f,'E:\Dropbox\Research\All Plottings\General plots\dynamicmodes_AVG-elastoviscousnumber_add20220216data.png','Resolution',500)

% buckling types vs. elastoviscous numner 
figure('color', 'w'); set(gcf, 'Position', [100 300 800 400]);
for no_Group = 1:NumGroup
    for buk = 1:2
        averaged_mu_bar_buk(no_Group, buk) = mean(DrawXY_buckletypes{no_Group, buk}(:, 2), 'omitnan');
    end
end
for PAtypeloop = 1:max(PAtype_label)  % Loop to draw according to the PA type labels.
    tmp5 = find(PAtype_label==PAtypeloop);
    Todraw_averaged_mu_bar_buk(PAtypeloop, :) = mean(averaged_mu_bar_buk(tmp5, :),1);
    h_buk_avg(PAtypeloop) = semilogx(Todraw_averaged_mu_bar_buk(PAtypeloop, :), 1:buk,'Color', CM(PAtypeloop,:), ...
        'LineStyle', 'none', 'marker', MyMarkers{1, PAtypeloop}, 'MarkerSize', 10,'LineWidth', 2); hold on
end
ylim([0 3]);
names = {'Global buckling'; 'Partial buckling'};
set(gca,'ytick',[1:2],'yticklabel',names,'FontSize',10);ytickangle(90)
xlabel('$\bar{\mu}$','FontSize', 15,'Interpreter', 'latex');
ylabel('$Dynamic\ Modes$','FontSize', 14,'Interpreter', 'latex');
legend(h_buk_avg, '$Square\ \lambda/a=1.88$' , '$Square\ \lambda/a=1.56$',...
    '$Rhombic\ \lambda/a=1.50$', '$Rhombic\ \lambda/a=1.75$','$Square\ \lambda/a=1.50$',...
    '$Square\ \lambda/a=1.50\ Flow\ angle=10^{\circ}$','Interpreter', 'latex', 'Location','best');  
%%% change the legends according to the excel 'ForActinPostprocessing.xlsx - sheet2'.
% f=gcf;
% exportgraphics(f,'E:\Dropbox\Research\All Plottings\General plots\bucklingtypes_AVG-elastoviscousnumber_add20220216data.png','Resolution',500)
