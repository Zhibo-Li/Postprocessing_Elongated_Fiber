% Postprocess code for the actin filaments in porous media (all in one).

%% Get the information, like contour length and elastoviscous number, from the reconstructed shape of the filaments.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
xlsfile = readcell('ForActinPostprocessing_laptop.xlsx','Sheet','Sheet1','NumHeaderLines',1); % This is the file that contains all the information about the later processing (in sheet 1).

NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.
PAsPath = xlsfile(:, 3);  % Path of the pillar array information.
savePath = xlsfile(:, 4);  % Saving path.
ChannelH = xlsfile(:, 5);  % Channel height (um)
ChannelW = xlsfile(:, 6);  % Channel width (um)
FlowRate = xlsfile(:, 7);  % Flow rate (nL/s)
Init_U = xlsfile(:, 8);  % Initial velocity (m/s)
Obj_Mag = xlsfile(:, 9); % Calibration (um/pixel)
C2C_dis = xlsfile(:, 12); % Center-to-center distance (um)
Pillar_a = xlsfile(:, 13); % Pillar diameter (um)

for no_Group = 1: NumGroup 
    
    Info(no_Group).ExpDate = ExpDate{no_Group, 1};  % list of the experiment date.
    filelist = dir(fullfile(storePath{no_Group},'*.mat'));  % list of the .mat files which contain the reconstruction information (came from 'Filaments detection' code) in one group.
    Info(no_Group).filelist = filelist;
    
    for no_Case = 1:length(filelist)

        load([storePath{no_Group}, filesep , filelist(no_Case).name])
%         lzero = max(lobject,ceil(5*lnoise));   % This is so important!!!!!!!!!!! Came from 'elongated_objects_in_flow.m'.
        lzero = 0;

        if Good_case(end) > xy.nframe  % This condition is not enough !!!
            Good_case_abs = Good_case;
            Good_case_frm = find(ismember(xy(1).frame, Good_case));
        else
            Good_case_frm = Good_case;
        end
        % depending on how to get the *.mat, Good_case could have different
        % meanings.
        
        sorted_lengths = sort(xy.arclen_spl(Good_case_frm));  % Good_case: index of the 'good' cases
%         contour_length = mean(sorted_lengths(round(numel(sorted_lengths)/10):end)) * Obj_Mag{no_Group};  % Select the 10% lengest filaments and averaged as the contour length. (UNIT: um)
        contour_length = max(sorted_lengths) * Obj_Mag{no_Group};  % Select the lengest filaments as the contour length (UNIT: um).
        Info(no_Group).contour_length(no_Case) = contour_length; 
        Info(no_Group).elastoviscousNum(no_Case) = VicFc_Get_elastoviscousNum(contour_length*1e-6, Init_U{no_Group}, Pillar_a{no_Group}*1e-6);  % Calculate the 'global' elastoviscous number for one trajectory.  (*1e-6: convert the unit to m)
        
        centroidxy = reshape(cell2mat(xy.centroid),2,numel(xy.centroid));
        centroidxy = centroidxy(:, Good_case_frm); 
%         centroidxy = centroidxy + lzero;  % The real CoM in the image!! BUT cannot add lzero here because following 'centroidxy_k' need it (NOT in IMAGE)!!!     
%         Amp(no_Case) = max(centroidxy(2,:)) -  min(centroidxy(2,:));  % This is to calculate the amplitude of the trajectory.


        centroidcom = complex(centroidxy(1,:), centroidxy(2,:)); % This is added on 03/12/2021 for calculating the minimum distance between CoM and pillar centres.
        centroidcom_shift = complex(lzero, lzero) + centroidcom;  % shift the coordinate  % Should notice here that need to shift both along x and y direction.
        StartPoint = centroidcom_shift(1);
        load(PAsPath{no_Group}); % Load the pillar array information.
        [~,idx] = sort(centers(:,1));
        sortedcenters = centers(idx,:);  % Sort the centres
        sortedcenters(:,2) = 2048 - sortedcenters(:,2); % Flip the coordinate because these are gotten from image coordinate system.
        sortedcenterscom = complex(sortedcenters(:,1), sortedcenters(:, 2));
        DistancesC2C_xy = StartPoint - sortedcenterscom; % The distances between the centers of the pillars and the COM of the filaments.
        [Ra,Ind] = min(DistancesC2C_xy);  % Ra: the minimum distance to one of the pillars (P), including x & y directions.     Ind: the index of the pillar (P).      Index of 'Ind': the index of the filaments.
        Info(no_Group).mini_Dis(no_Case) = abs(Ra) * Obj_Mag{no_Group};  % mini_Dis is the minimum distance between CoM and pillar centres (UNIT: um).
        

        centroidxy_plot_ind = 1;  % Plot -- x axis.
        count = 1;  % for plot
        for no_Goodcas = 1:size(Good_case,2)
            
            if no_Goodcas > 1
                if Good_case_frm(no_Goodcas) - TheLastOne >= 0  % Equal interval sampling. If put 4, it will be 5 and so on.      
                    centroidxy_plot_ind(end+1) = no_Goodcas;
                    
                    crd = xy.crd{1,Good_case_frm(no_Goodcas)};
                    centroidxy_k = centroidxy(:, no_Goodcas);
                    
                    L_ee = sqrt((crd(1,1)-crd(end,1))^2 + (crd(1,2)-crd(end,2))^2) * Obj_Mag{no_Group};  % Calculate the L_ee. (end-to-end distance)
                    Info(no_Group).L_ee_norm{no_Case}(count) = min(L_ee / contour_length, 1);  % Normalized by the contour length calculated above. 
                    % NOTE: if Lee/L0 is larger than 1 then Lee/L0 = 1. 
                    
                    Gyr = 1/size(crd,1) * [sum((crd(:, 1)-centroidxy_k(1)).^2),  sum((crd(:, 1)-centroidxy_k(1)) .* (crd(:, 2)-centroidxy_k(2)));
                        sum((crd(:, 2)-centroidxy_k(2)) .* (crd(:, 1)-centroidxy_k(1))), sum((crd(:, 2)-centroidxy_k(2)).^2)];
                    
                    [eigenV,eigenD] = eig(Gyr);
                    [d,ind] = sort(diag(eigenD));
                    Ds = eigenD(ind,ind);
                    Vs = eigenV(:,ind);
                    Info(no_Group).Chi{no_Case}(count) = atan(Vs(2,2)/Vs(1,2));
                    Info(no_Group).Chi_norm{no_Case}(count) = atan(Vs(2,2)/Vs(1,2))/pi;
                    
                    Lambda1 = eigenD(2,2); Lambda2 =  eigenD(1,1);
                    Info(no_Group).aniso{no_Case}(count) = 1 - 4*Lambda1*Lambda2/(Lambda1+Lambda2)^2;
                    
                    Info(no_Group).SelectCaseNo{no_Case}(count) = Good_case_frm(no_Goodcas);  % The indexes of the Goodcases. Mostly for the frequency calculation.
                    Info(no_Group).cont_lens_frm{no_Case}(count) = xy.arclen_spl(Good_case_frm(no_Goodcas)) * Obj_Mag{no_Group}; % The contour lengths in each frame (UNIT: um)!

                    count = count + 1;
                    TheLastOne = Good_case_frm(no_Goodcas);
                end
                
            else
                crd = xy.crd{1,Good_case_frm(no_Goodcas)};
                centroidxy_k = centroidxy(:, no_Goodcas);
                
                L_ee = sqrt((crd(1,1)-crd(end,1))^2 + (crd(1,2)-crd(end,2))^2) * Obj_Mag{no_Group};  % Calculate the L_ee. (end-to-end distance)
                Info(no_Group).L_ee_norm{no_Case}(count) = min(L_ee / contour_length, 1);  % Normalized by the contour length calculated above.  
                % NOTE: if Lee/L0 is larger than 1 then Lee/L0 = 1. 
                
                Gyr = 1/size(crd,1) * [sum((crd(:, 1)-centroidxy_k(1)).^2),  sum((crd(:, 1)-centroidxy_k(1)) .* (crd(:, 2)-centroidxy_k(2)));
                    sum((crd(:, 2)-centroidxy_k(2)) .* (crd(:, 1)-centroidxy_k(1))), sum((crd(:, 2)-centroidxy_k(2)).^2)];
                
                [eigenV,eigenD] = eig(Gyr);
                [d,ind] = sort(diag(eigenD));
                Ds = eigenD(ind,ind);
                Vs = eigenV(:,ind);
                Info(no_Group).Chi{no_Case}(no_Goodcas) = atan(Vs(2,2)/Vs(1,2));
                Info(no_Group).Chi_norm{no_Case}(count) = atan(Vs(2,2)/Vs(1,2))/pi;
                
                Lambda1 = eigenD(2,2); Lambda2 =  eigenD(1,1);
                Info(no_Group).aniso{no_Case}(no_Goodcas) = 1 - 4*Lambda1*Lambda2/(Lambda1+Lambda2)^2;
                
                Info(no_Group).SelectCaseNo{no_Case}(count) = Good_case_frm(no_Goodcas);  % The indexes of the Goodcases. Mostly for the frequency calculation.
                Info(no_Group).cont_lens_frm{no_Case}(count) = xy.arclen_spl(Good_case_frm(no_Goodcas)) * Obj_Mag{no_Group}; % The contour lengths in each frame (UNIT: um)!
                
                count = count + 1;
                TheLastOne = Good_case_frm(no_Goodcas);
            end
        end
        
        Info(no_Group).centroidxy_plotX{no_Case} = centroidxy(1,centroidxy_plot_ind) + lzero;  % lzzero is so important in image plotting!! Namely, when you are drawing based on the raw figure.
        Info(no_Group).centroidxy_plotY{no_Case} = centroidxy(2,centroidxy_plot_ind) + lzero;  % lzzero is so important!!
        
    end
    
end



%% Write the contour length and the elastoviscous numbers into the excel. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The excels are stored in ...Processing\EXP DATE\results normally.
% This part should run after the above 'Get the information' part.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for no_Group = 13: NumGroup

    disp(['You are processing ', datestr(Info(no_Group).ExpDate), ' data!']);
    [filename, pathname] = uigetfile('E:\Dropbox\Research\All Plottings\Manually classification documents\.xlsx', 'Please choose the excel accordingly and carefully!!');
    xlsfile = readcell([pathname, filename],'Sheet','Sheet2','NumHeaderLines',1);     % The files. (e.g. Classification manually 20220104-Actin.xlsx)
    thefiles = xlsfile(:, 14);

    filelist = dir(fullfile(storePath{no_Group},'*.mat'));  % list of the .mat files which contain the reconstruction information (came from 'Filaments detection' code) in one group.
    for no_Case = 1:length(filelist)
        name = Info(no_Group).filelist(no_Case).name(12:end-11);
        Ind = find(contains(thefiles,name));  % Search for text that has 'name' as part of the text.

        if isempty(Ind)
            disp('There are some problems with the alignment. Please check!!!');
        end

        Value = [Info(no_Group).contour_length(no_Case), Info(no_Group).elastoviscousNum(no_Case)];  % The contour length and the elastoviscous numbers.
        Loc = ['J', num2str(Ind+1), ':K', num2str(Ind+1)];  % The locations in the excel should be written into. (+1 because there is headerline in the excel.)
        writematrix(Value,[pathname, filename],'Sheet','Sheet2','Range', Loc);  % Write the value inti the excel.
    end
end



%% Plot 'L_ee_norm', 'omega' & 'Chi_norm'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part is to plot the characterizations along one trajectory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except Info xlsfile
NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.
PAsPath = xlsfile(:, 3);  % Path of the pillar array information.
savePath = xlsfile(:, 4);  % Saving path.
ChannelH = xlsfile(:, 5);  % Channel height (um)
ChannelW = xlsfile(:, 6);  % Channel width (um)
FlowRate = xlsfile(:, 7);  % Flow rate (nL/s)
Init_U = xlsfile(:, 8);  % Initial velocity (m/s)
Obj_Mag = xlsfile(:, 9); % Calibration (um/pixel)
C2C_dis = xlsfile(:, 12); % Center-to-center distance (um)
Pillar_a = xlsfile(:, 13); % Pillar diameter (um)

for no_Group = 1: NumGroup 
    
    filelist = dir(fullfile(storePath{no_Group},'*.mat'));  % list of the .mat files which contain the reconstruction information (came from 'Filaments detection' code) in one group.
    Info(no_Group).filelist = filelist;
    
    for no_Case = 1:length(filelist)
    figure('color', 'w'); set(gcf, 'Position', [100 300 1500 300]);
    
%     yyaxis left  % Plot 'omega'
%     plot(Info(no_Group).centroidxy_plotX{1,no_Case}, Info(no_Group).aniso{1,no_Case}, '--o', 'linewidth', 2);
%     xlim([0 2048]); ylim([0 1]);
%     xlabel('Position (px)','FontSize', 14,'Interpreter', 'latex')
%     ylabel('$\omega = 1 - 4{\lambda}_1{\lambda}_2/({\lambda}_1+{\lambda}_2)^2$','FontSize', 14,'Interpreter', 'latex');

%     hTxt(2,1)=text(-90,-0.3,{'$L_{ee}/L_0$'},'Interpreter', 'Latex','Rotation',0,'FontSize',20,'color','m'); 
%     hTxt(1,1)=text(-90,0.3,{'$\omega$'},'Interpreter', 'Latex','Rotation',0,'FontSize',20,'color','b');set(hTxt,'Rotation',90);
%     set(gca,'xtick',[])
%     Added @ 2021-11-24 to draw omega, Lee and Chi together!!
     
    yyaxis left  % Plot 'L_ee_norm'
    plot(Info(no_Group).centroidxy_plotX{1,no_Case}, Info(no_Group).L_ee_norm{1,no_Case}, '--o', 'linewidth', 2);
    xlim([0 2048]); ylim([0 max(Info(no_Group).L_ee_norm{1,no_Case})]);
    ylabel('$L_{ee}/L_0$','FontSize', 14,'Interpreter', 'latex');

    yyaxis right  % Plot 'Chi_norm'
    plot(Info(no_Group).centroidxy_plotX{1,no_Case}, Info(no_Group).Chi_norm{1,no_Case}, '--o', 'linewidth', 2);
    xlim([0 2048]); ylim([-0.5 0.5]);
    ylabel('${\chi} / {\pi}$','FontSize', 14,'Interpreter', 'latex');
    line([0,2048],[0,0],'LineStyle','--','Color',[0.8500 0.3250 0.0980],'LineWidth',0.2);  
    
    % This is to calculate the average positions of the pillars.
    [ave_PA_Ra, pillar_column, ~] = VicFc_Get_PAsInfo(PAsPath{no_Group});

    % The centers of the pillars
    for n_centers = 1:numel(pillar_column)
        if mod(n_centers, 2) == 0
            line([pillar_column(n_centers), pillar_column(n_centers)],[-0.5,0.5],...
                'Color','r','LineWidth',0.5);
        else
            line([pillar_column(n_centers), pillar_column(n_centers)],[-0.5,0.5],...
                'Color','c','LineWidth',0.5);  % Draw the positions of the circle's center
        end
    end
    
    % The edges of the pillars
    for n_centers = 1:numel(pillar_column)
        if mod(n_centers, 2) == 0
            line([pillar_column(n_centers)-ave_PA_Ra, pillar_column(n_centers)-ave_PA_Ra],...
                [-0.5,0.5],'LineStyle',':','Color','r','LineWidth',0.3);
            line([pillar_column(n_centers)+ave_PA_Ra, pillar_column(n_centers)+ave_PA_Ra],...
                [-0.5,0.5],'LineStyle',':','Color','r','LineWidth',0.3);
        else
            line([pillar_column(n_centers)-ave_PA_Ra, pillar_column(n_centers)-ave_PA_Ra],...
                [-0.5,0.5],'LineStyle',':','Color','c','LineWidth',0.3);
            line([pillar_column(n_centers)+ave_PA_Ra, pillar_column(n_centers)+ave_PA_Ra],...
                [-0.5,0.5],'LineStyle',':','Color','c','LineWidth',0.3);  % Draw the positions of the circle's edges.
        end
    end
    
%     export_fig([savePath{no_Group},filesep,Info(no_Group).filelist(no_Case).name,'_Lee_Chi_Position'],'-tif')
%     savefig([savePath{no_Group},filesep,Info(no_Group).filelist(no_Case).name,'_Lee_Chi_Position.fig'])

    end
end



%% Plot 'L_ee_norm', 'omega' & 'Chi_norm' in cell;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part is to plot the characterizations in one cell.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except Info xlsfile
NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.
PAsPath = xlsfile(:, 3);  % Path of the pillar array information.
savePath = xlsfile(:, 4);  % Saving path.
ChannelH = xlsfile(:, 5);  % Channel height (um)
ChannelW = xlsfile(:, 6);  % Channel width (um)
FlowRate = xlsfile(:, 7);  % Flow rate (nL/s)
Init_U = xlsfile(:, 8);  % Initial velocity (m/s)
Obj_Mag = xlsfile(:, 9); % Calibration (um/pixel)
C2C_dis = xlsfile(:, 12); % Center-to-center distance (um)
Pillar_a = xlsfile(:, 13); % Pillar diameter (um)

figure('color', 'w'); set(gcf, 'Position', [100 300 1500 300]);
for no_Group = 7:10%NumGroup
    
    % This is to calculate the average positions of the pillars.
    [~, pillar_column, ~] = VicFc_Get_PAsInfo(PAsPath{no_Group});
    DD = mean(diff(pillar_column));   % the distance betweeen two columns
    pillar_column_add = [pillar_column(1)-DD, pillar_column];
    
%     figure('color', 'w'); set(gcf, 'Position', [100 300 1500 300]);
    filelist = dir(fullfile(storePath{no_Group},'*.mat'));  % list of the .mat files which contain the reconstruction information (came from 'Filaments detection' code) in one group.
    Info(no_Group).filelist = filelist;
    
    for no_Case = 1:length(filelist)

        Cell_Position_Case = [];  % Cell_***_Case: superposition the *** in the cell for every case. (*** = Position, Lee, Chi). Here is to initialize for every case. ('Cell' here means one/half period)
        Cell_Lee_Case = [];
        Cell_Chi_Case = [];
        
%         figure('color', 'w'); set(gcf, 'Position', [-1600 300 1500 300]);
%         if Info(no_Group).contour_length(no_Case)>10; continue; end
%         if minRa(no_Case)<11.4050; continue; end    % The 'minRa' comes from the code 'Vic_Compare_Streamlines_to_Simulation.m'. To be continued...
        
%         for kk = 1:2:numel(pillar_column_add)-2  % every two columns
        for kk = 1:numel(pillar_column_add)-1  % every column
%             tmp_ind = Info(no_Group).centroidxy_plotX{1,no_Case} >= pillar_column_add(kk) & Info(no_Group).centroidxy_plotX{1,no_Case} <= pillar_column_add(kk+2);  % every two columns
            tmp_ind = Info(no_Group).centroidxy_plotX{1,no_Case} >= pillar_column_add(kk) & Info(no_Group).centroidxy_plotX{1,no_Case} <= pillar_column_add(kk+1);  % every column
            
%             yyaxis left
% %             plot((Info(no_Group).centroidxy_plotX{1,no_Case}(tmp_ind)-pillar_column_add(kk))/mean(diff(pillar_column_add))/2, Info(no_Group).aniso{1,no_Case}(tmp_ind),'marker','*','LineStyle', 'none'); hold on    % every two columns
%             plot((Info(no_Group).centroidxy_plotX{1,no_Case}(tmp_ind)-pillar_column_add(kk))/mean(diff(pillar_column_add)), Info(no_Group).aniso{1,no_Case}(tmp_ind),'marker','*','LineStyle', 'none'); hold on  % every column
%             xlabel('$Normalized\ Position\ during\ One\ Period$','FontSize', 14,'Interpreter', 'latex')
%             ylabel('$\omega = 1 - 4{\lambda}_1{\lambda}_2/({\lambda}_1+{\lambda}_2)^2$','FontSize', 14,'Interpreter', 'latex');
            
            yyaxis left
%             plot((Info(no_Group).centroidxy_plotX{1,no_Case}(tmp_ind)-pillar_column_add(kk))/mean(diff(pillar_column_add))/2, Info(no_Group).L_ee_norm{1,no_Case}(tmp_ind),'marker','*','LineStyle', 'none'); hold on    % every two columns
            plot((Info(no_Group).centroidxy_plotX{1,no_Case}(tmp_ind)-pillar_column_add(kk))/mean(diff(pillar_column_add)), Info(no_Group).L_ee_norm{1,no_Case}(tmp_ind),'marker','*','LineStyle', 'none'); hold on  % every column
            xlabel('$Normalized\ Position\ during\ One\ Period$','FontSize', 14,'Interpreter', 'latex') 
            ylabel('$L_{ee}/L_0$','FontSize', 14,'Interpreter', 'latex');
            
            yyaxis right
%             plot((Info(no_Group).centroidxy_plotX{1,no_Case}(tmp_ind)-pillar_column_add(kk))/mean(diff(pillar_column_add))/2, Info(no_Group).Chi_norm{1,no_Case}(tmp_ind),'marker','o','LineStyle', 'none');hold on    % every two columns
%             plot((Info(no_Group).centroidxy_plotX{1,no_Case}(tmp_ind)-pillar_column_add(kk))/mean(diff(pillar_column_add)), Info(no_Group).Chi_norm{1,no_Case}(tmp_ind),'marker','o','LineStyle', 'none');hold on  % every column
            plot((Info(no_Group).centroidxy_plotX{1,no_Case}(tmp_ind)-pillar_column_add(kk))/mean(diff(pillar_column_add)), abs(Info(no_Group).Chi_norm{1,no_Case}(tmp_ind)),'marker','o','LineStyle', 'none');hold on  % every column & absolute value
            xlim([0 1]); ylim([0 0.5]);
            ylabel('${\chi} / {\pi}$','FontSize', 14,'Interpreter', 'latex');

            % The following is to fit the data (Lee and Chi) in one/half
            % cell and check the reason for the 'streaks' occurring in the plot.
            Cell_Position_Case = [(Info(no_Group).centroidxy_plotX{1,no_Case}(tmp_ind)-pillar_column_add(kk))/mean(diff(pillar_column_add)), Cell_Position_Case];
            Cell_Lee_Case = [Info(no_Group).L_ee_norm{1,no_Case}(tmp_ind), Cell_Lee_Case];   
            Cell_Chi_Case = [abs(Info(no_Group).Chi_norm{1,no_Case}(tmp_ind)), Cell_Chi_Case];    

        end
%         FOO = [Cell_Position_Case; Cell_Lee_Case; Cell_Chi_Case];
%         Info(no_Group).Cell_Position_Lee_Chi_Case{1,no_Case} = sortrows(FOO')';  % Save the sorted 'position', 'Lee', 'Chi' for each case into Info.
% 
%         fit_gauss = fit(Info(no_Group).Cell_Position_Lee_Chi_Case{1,no_Case}(1,:).',Info(no_Group).Cell_Position_Lee_Chi_Case{1,no_Case}(3,:).','gauss1');

%         plot(fit_gauss,Info(no_Group).Cell_Position_Lee_Chi_Case{1,no_Case}(1,:),Info(no_Group).Cell_Position_Lee_Chi_Case{1,no_Case}(3,:)); %hold on; legend off;
%         set(gca,'FontSize',24);
%         xlabel('$Normalized\ Position\ during\ Half\ Period$','FontSize', 28,'Interpreter', 'latex');
%         ylabel('${\chi} / {\pi}$','FontSize', 28,'Interpreter', 'latex');
%         title('$Deformation\ and\ Orientation\ in\ Half\ Period\ (one\ case)$','FontSize', 32,'Interpreter', 'latex');  % Draw the fitting for one case.

%         Info(no_Group).fit_centroid(1,no_Case) = fit_gauss.b1;  % the centroid (location)
%         Info(no_Group).fit_amplitude(1,no_Case) = fit_gauss.a1;  %  the amplitude  
%         Info(no_Group).fit_shape(1,no_Case) = fit_gauss.c1;  % related to the peak width
%         Con = confint(fit_gauss,0.95); Con_b1 = Con(: ,2);
%         Info(no_Group).fit_centroid_confint(1,no_Case) = Con_b1(2)-Con_b1(1);


%         title('$Deformation\ and\ Orientation\ of\ the\ Filaments\ in\ One\ Period\ (all\ cases)$','FontSize', 18,'Interpreter', 'latex')
        
%         export_fig([savePath{no_Group},filesep,'Deformation_of_the_Filaments_every_column_',datestr(ExpDate{no_Group}),'-',num2str(no_Group)],'-tif')
%         savefig([savePath{no_Group},filesep,'Deformation_of_the_Filaments_every_column_',datestr(ExpDate{no_Group}),'-',num2str(no_Group),'.fig'
    end    
end
% clear FOO
f=gcf;
exportgraphics(f,'E:\Dropbox\Research\All Plottings\General plots\Deformation_of_the_Filaments_in_One_Period_20211029-20220104.png','Resolution',500)



%% Plot deviation angle vs. contour length or elastoviscous Number at different flow angle.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('E:\Dropbox\tmp\INfO.mat'); % This part needs the variable 'info'.
contourL = []; mu_bar = []; Slope = []; filelist = []; XY_plot_traj = []; counter = 1;

% Angle = 10° and Angle = 20°
for no_Group = 13: 18
    mu_bar = [mu_bar, Info(no_Group).elastoviscousNum];
    contourL = [contourL, Info(no_Group).contour_length];
    filelist = [filelist, Info(no_Group).filelist'];
    for no_Case = 1:length(Info(no_Group).contour_length)
        XXX = Info(no_Group).centroidxy_plotX{1, no_Case};
        YYY = Info(no_Group).centroidxy_plotY{1, no_Case};
        XY = [XXX; YYY];
        XY_plot_traj{counter} = XY;
        Slope(counter) = VicFc_Get_Slope(XY, 5, -0.1, 500, 500, 1400);
        counter = counter + 1;
    end
end

% Angle = 0°
for no_Group = 7: 8
    mu_bar = [mu_bar, Info(no_Group).elastoviscousNum];
    contourL = [contourL, Info(no_Group).contour_length];
    filelist = [filelist, Info(no_Group).filelist'];
    for no_Case = 1:length(Info(no_Group).contour_length)
        XXX = Info(no_Group).centroidxy_plotX{1, no_Case};
        YYY = Info(no_Group).centroidxy_plotY{1, no_Case};
        XY = [XXX; YYY];
        XY_plot_traj{counter} = XY;
        [p,S] =   polyfit(XY(1, :), XY(2, :),1);
        [y_fit,delta] = polyval(p,XY(1, :),S);
        if mean(delta) < 10
            Slope(counter) = p(1);
        else
            Slope(counter) = nan;
        end
        XY_plot_traj_calculatedpiece{counter} = nan;
        counter = counter + 1;
    end
end


Ang10_caseNum = length(Info(13).contour_length) + length(Info(14).contour_length) + length(Info(15).contour_length);
Ang20_caseNum = length(Info(16).contour_length) + length(Info(17).contour_length) + length(Info(18).contour_length);
Ang0_caseNum = length(Info(8).contour_length) + length(Info(9).contour_length);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Drawing: contour Length vs deviation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('color', 'w'); set(gcf, 'Position', [100 100 800 500]);

plot(contourL(Ang10_caseNum+Ang20_caseNum+1:end), abs(atand(Slope(Ang10_caseNum+Ang20_caseNum+1:end))-atand(0.0092)), ...
    'Color','m', 'LineStyle', 'none', 'marker', 'o', 'MarkerSize', 10,'LineWidth', 2); hold on;
Ang0 = [contourL(Ang10_caseNum+Ang20_caseNum+1:end);atand(Slope(Ang10_caseNum+Ang20_caseNum+1:end))-atand(0.0092)];
Ang0(:, isnan(Ang0(2, :)))=[];
pp3 = polyfit(Ang0(1, :),Ang0(2, :),1);
f3 = polyval(pp3,Ang0(1, :));
Ang0_fit = [Ang0(1, :);f3];
Ang0_fit = sortrows(Ang0_fit',1);
plot(Ang0_fit(:, 1),abs(Ang0_fit(:, 2)),':m', 'LineWidth',2); hold on;

plot(contourL(1:Ang10_caseNum), abs(atand(Slope(1:Ang10_caseNum))-10.65), ...
    'Color','r', 'LineStyle', 'none', 'marker', 'o', 'MarkerSize', 10,'LineWidth', 2); hold on;
Ang10 = [contourL(1:70);atand(Slope(1:70))-10.65]; % minus the slope of the pillar array (10.65)!
Ang10(:, isnan(Ang10(2, :)))=[];
pp1 = polyfit(Ang10(1, :),Ang10(2, :),1);
f1 = polyval(pp1,Ang10(1, :));
Ang10_fit = [Ang10(1, :);f1];
Ang10_fit = sortrows(Ang10_fit',1);
plot(Ang10_fit(:, 1),abs(Ang10_fit(:, 2)),':r','LineWidth',2); hold on;

plot(contourL(Ang10_caseNum+1:Ang10_caseNum+Ang20_caseNum), abs(atand(Slope(Ang10_caseNum+1:Ang10_caseNum+Ang20_caseNum))-19.4), ...
    'Color','c', 'LineStyle', 'none', 'marker', 'o', 'MarkerSize', 10,'LineWidth', 2); hold on;
Ang20 = [contourL(Ang10_caseNum+1:Ang10_caseNum+Ang20_caseNum);atand(Slope(Ang10_caseNum+1:Ang10_caseNum+Ang20_caseNum))-19.4]; % minus the slope of the pillar array (19.4)!
Ang20(:, isnan(Ang20(2, :)))=[];
pp2 = polyfit(Ang20(1, :),Ang20(2, :),1);
f2 = polyval(pp2,Ang20(1, :));
Ang20_fit = [Ang20(1, :);f2];
Ang20_fit = sortrows(Ang20_fit',1);
plot(Ang20_fit(:, 1),abs(Ang20_fit(:, 2)),':c', 'LineWidth',2); hold on;

% title('CrMask Square Phi20 Gap10 FlowAng10', FontSize=12);
set(gca,'FontSize',16);
xlabel('$\bf{Contour\ length\ ({\mu}m)}$','FontSize', 20,'Interpreter', 'latex');
ylabel('$\bf{CoM\ trajectory\ deviation\ \phi\ (^{\circ})}$','FontSize', 20,'Interpreter', 'latex');
% hold on;plot([0 60], [0.1763 0.1763], 'b--', 'LineWidth',1.5)
legend({'Flow angle = 0°', '', 'Flow angle = 10°','', 'Flow angle = 20°',''},'FontSize', 20, 'Location','best')
ylim([0 10])
% f=gcf;
% exportgraphics(f,'E:\Dropbox\Research\All Plottings\General plots\20210914_20220216-20220217_trajectoryslope_contourlength_abs_degree.png','Resolution',100)

% 10 degree trajectory from simulation
data = readmatrix('D:\Dropbox\PROCESS remotely\202208_differentFlowangles_relatedto_0811exp_45deg\10deg\Data\Streamline.csv');
XXYY_Simu = data(1:10:end, 14:15); 
XXYY_Simu(XXYY_Simu(:, 2) > -0.00105, :) = []; XXYY_Simu(XXYY_Simu(:, 2) < -0.00140, :) = []; 
Slope10deg = VicFc_Get_Slope(sortrows(XXYY_Simu)', 5, 0, 0.0011, 0.0015, 0.0015);
hold on; plot([0 100], [10-atand(Slope10deg) 10-atand(Slope10deg)]);

data = readmatrix('D:\Dropbox\PROCESS remotely\202208_differentFlowangles_relatedto_0811exp_45deg\20deg\Data\Streamline.csv');
XXYY_Simu = data(1:10:end, 14:15); 
XXYY_Simu(XXYY_Simu(:, 2) > 0, :) = []; XXYY_Simu(XXYY_Simu(:, 2) < -0.0006, :) = []; 
Slope20deg = VicFc_Get_Slope(sortrows(XXYY_Simu)', 5, 0.01, 0.0005, 0.002, 0.002);
hold on; plot([0 100], [20-atand(Slope20deg) 20-atand(Slope20deg)])



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Drawing: elastoviscous Number vs deviation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('color', 'w'); set(gcf, 'Position', [100 100 800 500]);

semilogx(mu_bar(Ang10_caseNum+Ang20_caseNum+1:end), abs(atand(Slope(Ang10_caseNum+Ang20_caseNum+1:end))-atand(0.0092)), ...
    'Color','m', 'LineStyle', 'none', 'marker', 'o', 'MarkerSize', 10,'LineWidth', 2); hold on;
Ang0 = [log10(mu_bar(Ang10_caseNum+Ang20_caseNum+1:end));atand(Slope(Ang10_caseNum+Ang20_caseNum+1:end))-atand(0.0092)];
Ang0(:, isnan(Ang0(2, :)))=[];
pp3 = polyfit(Ang0(1, :),Ang0(2, :),1);
f3 = polyval(pp3,Ang0(1, :));
Ang0_fit = [Ang0(1, :);f3];
Ang0_fit = sortrows(Ang0_fit',1);
plot(10.^(Ang0_fit(:, 1)),abs(Ang0_fit(:, 2)),':m', 'LineWidth',2); hold on;

semilogx(mu_bar(1:Ang10_caseNum), abs(atand(Slope(1:Ang10_caseNum))-10), ...
    'Color','r', 'LineStyle', 'none', 'marker', 'o', 'MarkerSize', 10,'LineWidth', 2); hold on;
Ang10 = [log10(mu_bar(1:70));atand(Slope(1:70))-10]; % minus the slope of the pillar array!
Ang10(:, isnan(Ang10(2, :)))=[];
pp1 = polyfit(Ang10(1, :),Ang10(2, :),1);
f1 = polyval(pp1,Ang10(1, :));
Ang10_fit = [Ang10(1, :);f1];
Ang10_fit = sortrows(Ang10_fit',1);
plot(10.^(Ang10_fit(:, 1)),abs(Ang10_fit(:, 2)),':r','LineWidth',2); hold on;

semilogx(mu_bar(Ang10_caseNum+1:Ang10_caseNum+Ang20_caseNum), abs(atand(Slope(Ang10_caseNum+1:Ang10_caseNum+Ang20_caseNum))-20), ...
    'Color','c', 'LineStyle', 'none', 'marker', 'o', 'MarkerSize', 10,'LineWidth', 2); hold on;
Ang20 = [log10(mu_bar(Ang10_caseNum+1:Ang10_caseNum+Ang20_caseNum));atand(Slope(Ang10_caseNum+1:Ang10_caseNum+Ang20_caseNum))-20];
Ang20(:, isnan(Ang20(2, :)))=[];
pp2 = polyfit(Ang20(1, :),Ang20(2, :),1);
f2 = polyval(pp2,Ang20(1, :));
Ang20_fit = [Ang20(1, :);f2];
Ang20_fit = sortrows(Ang20_fit',1);
plot(10.^(Ang20_fit(:, 1)),abs(Ang20_fit(:, 2)),':c', 'LineWidth',2); hold on;

% title('CrMask Square Phi20 Gap10 FlowAng10', FontSize=12);
set(gca,'FontSize',16);
xlabel('$\bf{Elastoviscous\ Number\ \bar{\mu}}$','FontSize', 20,'Interpreter', 'latex');
ylabel('$\bf{CoM\ trajectory\ deviation\ \phi\ (^{\circ})}$','FontSize', 20,'Interpreter', 'latex');
% hold on;plot([0 60], [0.1763 0.1763], 'b--', 'LineWidth',1.5)
legend({'Flow angle = 0°', '', 'Flow angle = 10°','', 'Flow angle = 20°',''},'FontSize', 20, 'Location','best')
ylim([0 10])
% f=gcf;
% exportgraphics(f,'E:\Dropbox\Research\All Plottings\General plots\20210914_20220216-20220217_trajectoryslope_elastoviscousnumber_abs_degree.png','Resolution',100)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Pick a point on 'elastoviscous Number vs deviation' and draw the trajectory as well as find the video of the case.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[the_elsvisN, the_devi] = ginput(1); % pick up the point you want to show the trajectory.
% [the_conL, the_devi] = ginput(1);

the_layout = 10; % declare the layout (flow angle) for the converting later.

% the_loc = intersect(find(contourL>the_conL*0.99 & contourL<the_conL*1.01),...
%     find(abs(atand(Slope)-the_layout)>the_devi*0.99 & abs(atand(Slope)-the_layout)...
%     < the_devi*1.01));

the_loc = intersect(find(mu_bar>the_elsvisN*0.95 & mu_bar<the_elsvisN*1.05),...
    find(abs(atand(Slope)-the_layout)>the_devi*0.95 & abs(atand(Slope)-the_layout)...
    < the_devi*1.05));  % Change the control range if there is an error.

filelist(the_loc).name     % show the file
filelist(the_loc).folder    % show the folder

index_tmp = find(strcmp(xlsfile(:, 2), filelist(the_loc).folder)); % identify the PA information (Path of the pillar array).
PAsPath = xlsfile(index_tmp(1), 3);

figure('color', 'w', 'Position',[100 100 800 600]); 
load(PAsPath{1}); % load the PAs information
centers = cat(2, centers(:,1), 2048 - centers(:,2)); % 2048 is image size.
viscircles(centers,radii,'LineStyle','--', 'LineWidth', 1, 'Color', 'k'); hold on

plot(XY_plot_traj{the_loc}(1, :), XY_plot_traj{the_loc}(2, :), 'LineWidth',2); axis equal; hold on
plot(XY_plot_traj_calculatedpiece{the_loc}(1, :), XY_plot_traj_calculatedpiece{the_loc}(2, :),'r','LineWidth',2); hold off
axis off


%% Plot elastoviscousNum vs. L_ee_morm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part is used to polt the elasto-viscous number /mu vs. end-to-end
% length L_ee_norm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except Info xlsfile
NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.
PAsPath = xlsfile(:, 3);  % Path of the pillar array information.
savePath = xlsfile(:, 4);  % Saving path.
ChannelH = xlsfile(:, 5);  % Channel height (um)
ChannelW = xlsfile(:, 6);  % Channel width (um)
FlowRate = xlsfile(:, 7);  % Flow rate (nL/s)
Init_U = xlsfile(:, 8);  % Initial velocity (m/s)
Obj_Mag = xlsfile(:, 9); % Calibration (um/pixel)
C2C_dis = xlsfile(:, 12); % Center-to-center distance (um)
Pillar_a = xlsfile(:, 13); % Pillar diameter (um)

plot_counter1 = 1; counter1 = 1;  % plot_counter1: for the final plot. 
for no_Group = 1: NumGroup  
    
    for no_Case = 1:length(Info(no_Group).filelist)
        
        tmp_L_ee_norm = Info(no_Group).L_ee_norm{1, no_Case};  % The L_ee_norm for one trajectory.
         
        [values,BinEdges] = histcounts(tmp_L_ee_norm, 10);  % Histogram bin counts. The number of bins is 10.
        tmp_ind = find(values==max(values));  % The maximum after histogram.
        if numel(tmp_ind) > 1  % The cases which don't only have one maximum value.
            MultiMaxCases{counter1} = Info(no_Group).filelist(no_Case).name;
            counter1 = counter1 + 1;
        end
        L_ee_plot(plot_counter1) = (BinEdges(tmp_ind(1)+1) + BinEdges(tmp_ind(1)))/2;  % Most probable L_ee_norm.
%         L_ee_plot(plot_counter1) = (BinEdges(1) + BinEdges(2))/2;  % Minimum L_ee_norm.
        mu_0_plot(plot_counter1) = Info(no_Group).elastoviscousNum(no_Case);  % The elasto viscous number.
        plot_counter1 = plot_counter1 + 1;
        
    end
end

figure('color', 'w'); set(gcf, 'Position', [100 300 1000 500]);
semilogx(mu_0_plot, L_ee_plot, 'b*')
% title('$CoM\ Distribution$','FontSize', 18,'Interpreter', 'latex')
xlabel('$\mu_0$','FontSize', 18,'Interpreter', 'latex');
ylabel('$L_{ee}/L_0$','FontSize', 18,'Interpreter', 'latex');
% axis equal;
% xlim([0 2]);  ylim([0 1]);
% export_fig([savePath{no_Group},filesep,'mostprob_Lee_vs_mu0_',datestr(ExpDate{no_Group}),'-',num2str(no_Group)],'-tif')
% savefig([savePath{no_Group},filesep,'mostprob_Lee_vs_mu0_',datestr(ExpDate{no_Group}),'-',num2str(no_Group),'.fig'])



%%  Plot elastoviscousNum vs. L_ee_morm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw Lee/L0 vs mu0, meanwhile divide into different parts according to
% the initial position.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Plot pdf of CoM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part is used to polt the pdf of CoM in a cell.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except Info xlsfile
NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.
PAsPath = xlsfile(:, 3);  % Path of the pillar array information.
savePath = xlsfile(:, 4);  % Saving path.
ChannelH = xlsfile(:, 5);  % Channel height (um)
ChannelW = xlsfile(:, 6);  % Channel width (um)
FlowRate = xlsfile(:, 7);  % Flow rate (nL/s)
Init_U = xlsfile(:, 8);  % Initial velocity (m/s)
Obj_Mag = xlsfile(:, 9); % Calibration (um/pixel)
C2C_dis = xlsfile(:, 12); % Center-to-center distance (um)
Pillar_a = xlsfile(:, 13); % Pillar diameter (um)

for no_Group = 1: NumGroup
    
    filelist = dir(fullfile(storePath{no_Group},'*.mat'));  % list of the .mat files which contain the reconstruction information (came from 'Filaments detection' code) in one group.
    % This is to calculate the average positions of the pillars.
    [ave_PA_Ra, pillar_column, pillar_row] = VicFc_Get_PAsInfo(PAsPath{no_Group});
    DD = mean(diff(pillar_column));   % the distance betweeen two columns
    pillar_column_add = [pillar_column(1)-DD, pillar_column];
    Y_shift =  2048-pillar_row(end);   % shift the plot up or down to make it in a 'cell'. Need to be flipped and 2048 is the size of the image.
    
    for no_Case = 1:length(filelist)    % choose the cases to draw
        for kk = 1:2:numel(pillar_column_add)-2   % every two columns. 
            tmp_ind = Info(no_Group).centroidxy_plotX{1,no_Case} >= pillar_column_add(kk) & Info(no_Group).centroidxy_plotX{1,no_Case} <= pillar_column_add(kk+2);  % every two columns.
            scatter1 = scatter((Info(no_Group).centroidxy_plotX{1,no_Case}(tmp_ind)-pillar_column_add(kk))/DD, (mod(Info(no_Group).centroidxy_plotY{1,no_Case}(tmp_ind)-Y_shift,DD))/DD,'filled', 'b');
            alpha(scatter1,0.2);  hold on
            hold on;
        end
    end
    
    clearvars pillar_column_add pillar_column pillar_row
    
end

figure('color', 'w'); %set(gcf, 'Position', [100 300 1000 500]);
title('$CoM\ Distribution$','FontSize', 18,'Interpreter', 'latex')
viscircles([1,0; 1,1; 0,0.5; 2,0.5],[ave_PA_Ra/DD; ave_PA_Ra/DD; ave_PA_Ra/DD; ave_PA_Ra/DD]);  % plot the pillar array according to the layout.
xlabel('$x/D_{C2C}$','FontSize', 18,'Interpreter', 'latex');
ylabel('$y/D_{C2C}$','FontSize', 18,'Interpreter', 'latex');
axis equal;
xlim([0 2]);  ylim([0 1]);
% export_fig(['F:\PhD, PMMH, ESPCI\Processing\20210430-Actin\results\Figures\CoM_1nL'],'-tif')
% export_fig([savePath{no_Group},filesep,'CoM_',datestr(ExpDate{no_Group}),'-',num2str(no_Group)],'-tif')
% savefig([savePath{no_Group},filesep,'CoM_',datestr(ExpDate{no_Group}),'-',num2str(no_Group),'.fig'])



%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw Lee/L0 vs mu0, meanwhile divide into different parts according to
% the initial position.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% About the periodicity (Orientation & Deformation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except Info xlsfile
NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.
PAsPath = xlsfile(:, 3);  % Path of the pillar array information.
savePath = xlsfile(:, 4);  % Saving path.
ChannelH = xlsfile(:, 5);  % Channel height (um)
ChannelW = xlsfile(:, 6);  % Channel width (um)
FlowRate = xlsfile(:, 7);  % Flow rate (nL/s)
Init_U = xlsfile(:, 8);  % Initial velocity (m/s)
Obj_Mag = xlsfile(:, 9); % Calibration (um/pixel)
C2C_dis = xlsfile(:, 12); % Center-to-center distance (um)
Pillar_a = xlsfile(:, 13); % Pillar diameter (um)

plot_counter1 = 1;
for no_Group = 7:10%1: NumGroup

    The_GAP = C2C_dis{no_Group}*1e-6;  % *1e-6: convert the unit to m.
    filelist = dir(fullfile(storePath{no_Group},'*.mat'));  % list of the .mat files which contain the reconstruction information (came from 'Filaments detection' code) in one group.
    for no_Case = 1:length(filelist)    % choose the cases to draw
        [pxx,f] = plomb(Info(no_Group).Chi_norm{1, no_Case}, Info(no_Group).SelectCaseNo{1, no_Case},'power'); % Chi_norm
        [pxx1,f1] = plomb(Info(no_Group).L_ee_norm{1, no_Case}, Info(no_Group).SelectCaseNo{1, no_Case},'power'); % L_ee_norm
%         figure; plomb(Info(no_Group).L_ee_norm{1, no_Case}, Info(no_Group).SelectCaseNo{1, no_Case},'power');
        TMP = diff(Info(no_Group).Chi_norm{1, no_Case});
%         if min(cos(TMP)) < 0.9  % Select the smooth cases.
        if size(Info(no_Group).SelectCaseNo{1, no_Case},2) / max(Info(no_Group).SelectCaseNo{1, no_Case}) < 0.8  
            %         figure; plot(Info(no_Group).Chi_norm{1, no_Case});
            %         figure; plot(cos(TMP)); ylim([0.5 1]);
            Period_plot(plot_counter1) = 0;
            Period_plot1(plot_counter1) = 0;
        else
            %         figure; plot(pxx, f);
            [pk,f0] = findpeaks(pxx,f); % Chi_norm
            [pk1,f01] = findpeaks(pxx1,f1); % L_ee_norm
            if  isempty(f0)
                Period_plot(plot_counter1) = 0;  % Chi_norm
                Period_plot1(plot_counter1) = 0;   % L_ee_norm
            else
% % %                 figure('color', 'w'); set(gcf, 'Position', [100 300 1000 500]);
% % %                 plot(Info(no_Group).SelectCaseNo{1, no_Case}, Info(no_Group).Chi_norm{1, no_Case});set(gca,'FontSize',12)
% % %                 xlabel('$Time\ (frame)$','FontSize', 18,'Interpreter', 'latex');
% % %                 ylabel('${\chi} / {\pi}$','FontSize', 18,'Interpreter', 'latex');
% % % %                 export_fig([savePath{no_Group},filesep,'Orientation_vs_Position_',datestr(ExpDate{no_Group}),'-',num2str(no_Group)],'-tif')
% % %                 
% % %                 figure('color', 'w'); set(gcf, 'Position', [100 300 600 500]);
% % %                 plomb(Info(no_Group).Chi_norm{1, no_Case}, Info(no_Group).SelectCaseNo{1, no_Case},'power');set(gca,'FontSize',12)
% % % %                 export_fig([savePath{no_Group},filesep,'PSD_vs_Frequency_',datestr(ExpDate{no_Group}),'-',num2str(no_Group)],'-tif')
% % %                 
% % %                 figure('color', 'w'); set(gcf, 'Position', [100 300 600 500]);
% % %                 plot(f, pxx, f0, pk,'ro'); set(gca,'FontSize',12)
% % %                 xlabel('$Frequency\ (Hz)$','FontSize', 18,'Interpreter', 'latex');
% % %                 ylabel('$Power\ Spectral\ Density\ (PSD)$','FontSize', 18,'Interpreter', 'latex');
%                 export_fig([savePath{no_Group},filesep,'PSD_vs_Frequency2_',datestr(ExpDate{no_Group}),'-',num2str(no_Group)],'-tif')

                Period_plot(plot_counter1) = 1/f0(pk == max(pk)) * 0.02 * Init_U{no_Group, 1} / The_GAP; % Chi_norm.    Normalized  by the initial velocity and the c2c distanse between two pillars.
                Period_plot1(plot_counter1) = 1/f01(pk1 == max(pk1)) * 0.02 * Init_U{no_Group, 1} / The_GAP; % L_ee_norm
%                 if 1/f0(pk == max(pk)) > 90
%                     disp(Info(no_Group).filelist(no_Case).name)
%                 end
            end
        end
        mu_0_plot(plot_counter1) = Info(no_Group).elastoviscousNum(no_Case);  % The elasto viscous number.
        plot_counter1 = plot_counter1 + 1;
    end
    
end


PLOT = [mu_0_plot; Period_plot; Period_plot1];
PLOT(:,PLOT(2,:)==0)=[];
figure('color', 'w'); set(gcf, 'Position', [100 300 1000 500]);
semilogx(PLOT(1,:), PLOT(2,:), 'b*', 'MarkerSize', 10); hold on
semilogx(PLOT(1,:), PLOT(3,:), 'ro', 'MarkerSize', 10); hold on
% legend('$\chi$','FontSize', 18,'Interpreter', 'latex');


% semilogx(mu_0_plot, Period_plot, 'b*', 'MarkerSize', 10); hold on
% semilogx(mu_0_plot, Period_plot1, 'ro', 'MarkerSize', 10); hold on
set(gca,'FontSize',12)
% title('$CoM\ Distribution$','FontSize', 18,'Interpreter', 'latex')
xlabel('$\mu_0$','FontSize', 18,'Interpreter', 'latex');
ylabel('$U_0/(\lambda \times f_p)$','FontSize', 18,'Interpreter', 'latex');
% axis equal;
% xlim([0 2]);  
line([100 1e6], [1 1],'Color','red','LineStyle','--')
% legend('$\chi$','','FontSize', 18,'Interpreter', 'latex');
% ylim([0.4 1.4]);
% export_fig([savePath{no_Group},filesep,'Orientation_periodicity_vs_mu0_',datestr(ExpDate{no_Group}),'-',num2str(no_Group)],'-tif')
% savefig([savePath{no_Group},filesep,'mostprob_Lee_vs_mu0_',datestr(ExpDate{no_Group}),'-',num2str(no_Group),'.fig']) 
    
