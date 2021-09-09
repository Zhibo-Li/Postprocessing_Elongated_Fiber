clear; clc; close all;
xlsfile = readcell('ForActinPostprocessing.xlsx','NumHeaderLines',1); % This is the file contains all the information about the later processing.

NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
storePath = xlsfile(:, 2);  % Path of the data to be processed.
PAsPath = xlsfile(:, 3);  % Path of the pillar array information.
savePath = xlsfile(:, 4);  % Saving path.
ChannelH = xlsfile(:, 5);  % Channel height (um)
ChannelW = xlsfile(:, 6);  % Channel width (um)
FlowRate = xlsfile(:, 7);  % Flow rate (nL/s)
Init_U = xlsfile(:, 8);  % Initial velocity (m/s)
Obj_Mag = xlsfile(:, 9); % Calibration (um/pixel)
 

for no_Group = 1: NumGroup 
    
    filelist = dir(fullfile(storePath{no_Group},'*.mat'));  % list of the .mat files which contain the reconstruction information (came from 'Filaments detection' code) in one group.
    Info(no_Group).filelist = filelist;
    
    for no_Case = 1:length(filelist)
        
        load([storePath{no_Group}, filesep , filelist(no_Case).name])
        lzero = max(lobject,ceil(5*lnoise));   % This is so important!!!!!!!!!!! Came from when we do the filaments detection.
        
        sorted_lengths = sort(xy.arclen_spl(Good_case));  % Good_case: index of the 'good' cases
        contour_length = mean(sorted_lengths(round(numel(sorted_lengths)/10):end)) * Obj_Mag{no_Group};  % Select the 10% lengest filaments and averaged as the contour length. (UNIT: um)
        Info(no_Group).contour_length(no_Case) = contour_length; 
        Info(no_Group).elastoviscousNum(no_Case) = Get_elastoviscousNum(contour_length*1e-6, Init_U{no_Group});  % Calculate the 'global' elastoviscous number for one trajectory.  (*1e-6: convert the unit to m)
        
        centroidxy = reshape(cell2mat(xy.centroid),2,numel(xy.centroid));
        centroidxy = centroidxy(:, Good_case); 
        centroidxy = centroidxy + lzero;  % The CoM.
%         Amp(no_Case) = max(centroidxy(2,:)) -  min(centroidxy(2,:));  % This is to calculate the amplitude of the trajectory.
        
        centroidxy_plot_ind = 1;  % Plot -- x axis.
        count = 1;  % for plot
        for no_Goodcas = 1:size(Good_case,2)
            
            if no_Goodcas > 1
                if Good_case(no_Goodcas) - TheLastOne >= 4  % Equal interval sampling. If put 4, it will be 5 and so on.      
                    centroidxy_plot_ind(end+1) = no_Goodcas;
                    
                    crd = xy.crd{1,Good_case(no_Goodcas)};
                    centroidxy_k = centroidxy(:, no_Goodcas);
                    
                    L_ee = sqrt((crd(1,1)-crd(end,1))^2 + (crd(1,2)-crd(end,2))^2) * Obj_Mag{no_Group};  % Calculate the L_ee. (end-to-end distance)
                    Info(no_Group).L_ee_norm{no_Case}(count) = L_ee / contour_length;  % Normalized by the contour length calculated above.
                    
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
                    
                    count = count + 1;
                    TheLastOne = Good_case(no_Goodcas);
                end
                
            else
                crd = xy.crd{1,Good_case(no_Goodcas)};
                centroidxy_k = centroidxy(:, no_Goodcas);
                
                L_ee = sqrt((crd(1,1)-crd(end,1))^2 + (crd(1,2)-crd(end,2))^2) * Obj_Mag{no_Group};  % Calculate the L_ee. (end-to-end distance)
                Info(no_Group).L_ee_norm{no_Case}(count) = L_ee / contour_length;  % Normalized by the contour length calculated above.
                
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
                
                count = count + 1;
                
                TheLastOne = Good_case(no_Goodcas);
            end
        end
        
        Info(no_Group).centroidxy_plot{no_Case} = centroidxy(1,centroidxy_plot_ind);
        
    end
    
end



Info(1).filelist.name



%% Plot elastoviscousNum vs. L_ee_morm
% This part is used to polt the elasto-viscous number /mu vs. end-to-end
% length L_ee_norm.

clearvars -except Info

figure('color', 'w'); set(gcf, 'Position', [100 300 1000 500]);

fooo = 1; ffff = 1;
for no_Group = 1: NumGroup  
    
    for no_Case = 1:length(Info(no_Group).filelist)
        
        tmp_L_ee_norm = Info(no_Group).L_ee_norm{1, no_Group};
        
%         h = histogram(tmp_L_ee, 10);  
%         values = h.Values;   BinEdges = h.BinEdges;
        
        [values,BinEdges] = histcounts(tmp_L_ee_norm, 10);
        tmp_ind = find(values==max(values));
        if numel(tmp_ind) > 1
            yyyyy{ffff} = [num2str(no_Group), num2str(no_Group)];
            ffff = ffff + 1;
        end
        L_ee_plot(fooo) = (BinEdges(tmp_ind(1)+1) + BinEdges(tmp_ind(1)))/2;
%         L_ee_plot(fooo) = (BinEdges(1) + BinEdges(2))/2;
        mu_0_plot(fooo) = Info(no_Group).elastoviscousNum(no_Group);
        fooo = fooo + 1;
        
    end
end


semilogx(mu_0_plot, L_ee_plot, 'b*')
% scatter(mu_0_plot, L_ee_plot);
% title('$CoM\ Distribution$','FontSize', 18,'Interpreter', 'latex')

xlabel('$\mu_0$','FontSize', 18,'Interpreter', 'latex');
ylabel('$L_{ee}/L_0$','FontSize', 18,'Interpreter', 'latex');
% axis equal;
% xlim([0 2]);  ylim([0 1]);
% export_fig(['F:\PhD, PMMH, ESPCI\Processing\20210430-Actin\results\Figures\mostprob_Lee_vs_mu0'],'-tif')




for no_Group = 1: NumGroup 
    
    if no_Group == 1
        
        the_chosen_folder = '20210413-Actin';
        filelist=dir(fullfile('F:\PhD, PMMH, ESPCI\Processing',the_chosen_folder,'results','*.mat'));
        load('F:\PhD, PMMH, ESPCI\Processing\20210413-Actin\circlesforPAs1_S10.mat')
        ave_R = mean(radii);
        sorted_centers = sort(centers(:,1));
        Jumps = find(diff(sorted_centers) > 100); Jumps = [0;Jumps;size(sorted_centers,1)];
        for ii = 1:size(Jumps, 1)-1
            pillar_column(ii) = mean(sorted_centers(Jumps(ii)+1: Jumps(ii + 1)));
        end
        DD = mean(diff(pillar_column));   % the distance betweeen two columns
        pillar_column = [pillar_column(1)-DD, pillar_column];
        
%         sorted_centersY = sort(centers(:,2));
%         JumpsY = find(diff(sorted_centersY) > 80); JumpsY = [0;JumpsY;size(sorted_centersY,1)];
%         for ii = 1:size(JumpsY, 1)-1
%             pillar_row(ii) = mean(sorted_centers(JumpsY(ii)+1: JumpsY(ii + 1)));
%         end
%         Y_shift = pillar_row(1)/2;   % shift the plot up or down to make it in a 'cell'
        
        for jj = 1:length(filelist)    % choose the cases to draw
            for kk = 1:2:numel(pillar_column)-2  % in one period
                tmp_ind = Info(no_Group).centroidxy_plotX{1,jj} >= pillar_column(kk) & Info(no_Group).centroidxy_plotX{1,jj} <= pillar_column(kk+2);  % in one period
                %         plot((centroidxy_plotX{1,jj}(tmp_ind)-pillar_column(kk))/DD, (mod(centroidxy_plotY{1,jj}(tmp_ind)-Y_shift,DD))/DD,'Color','b','marker','.','LineStyle', 'none'); hold on    % in one period
                scatter1 = scatter((Info(no_Group).centroidxy_plotX{1,jj}(tmp_ind)-pillar_column(kk))/DD, (mod(Info(no_Group).centroidxy_plotY{1,jj}(tmp_ind)-Y_shift,DD))/DD,'filled', 'b');
                alpha(scatter1,0.2);  hold on
                hold on;
            end
        end
        
    else
        
        clearvars pillar_column pillar_row
        the_chosen_folder = '20210430-Actin';
        filelist=dir(fullfile('F:\PhD, PMMH, ESPCI\Processing',the_chosen_folder,'results','*.mat'));
        load('F:\PhD, PMMH, ESPCI\Processing\20210430-Actin\circlesforPAs2_S10.mat')
        ave_R = mean(radii);
        sorted_centers = sort(centers(:,1));
        Jumps = find(diff(sorted_centers) > 100); Jumps = [0;Jumps;size(sorted_centers,1)];
        for ii = 1:size(Jumps, 1)-1
            pillar_column(ii) = mean(sorted_centers(Jumps(ii)+1: Jumps(ii + 1)));
        end
        DD = mean(diff(pillar_column));   % the distance betweeen two columns
        pillar_column = [pillar_column(1)-DD, pillar_column];
        
%         sorted_centersY = sort(centers(:,2));
%         JumpsY = find(diff(sorted_centersY) > 80); JumpsY = [0;JumpsY;size(sorted_centersY,1)];
%         for ii = 1:size(JumpsY, 1)-1
%             pillar_row(ii) = mean(sorted_centers(JumpsY(ii)+1: JumpsY(ii + 1)));
%         end
%         Y_shift = pillar_row(1)/2;   % shift the plot up or down to make it in a 'cell'

        for jj = 1:21   % choose the cases to draw
            for kk = 1:2:numel(pillar_column)-2  % in one period
                tmp_ind = Info(no_Group).centroidxy_plotX{1,jj} >= pillar_column(kk) & Info(no_Group).centroidxy_plotX{1,jj} <= pillar_column(kk+2);  % in one period
                %         plot((centroidxy_plotX{1,jj}(tmp_ind)-pillar_column(kk))/DD, (mod(centroidxy_plotY{1,jj}(tmp_ind)-Y_shift,DD))/DD,'Color','b','marker','.','LineStyle', 'none'); hold on    % in one period
                scatter1 = scatter((Info(no_Group).centroidxy_plotX{1,jj}(tmp_ind)-pillar_column(kk))/DD, (mod(Info(no_Group).centroidxy_plotY{1,jj}(tmp_ind)-Y_shift,DD))/DD,'filled', 'b');
                alpha(scatter1,0.2);  hold on
                hold on;
            end
        end
    end
end

title('$CoM\ Distribution$','FontSize', 18,'Interpreter', 'latex')
viscircles([1,0; 1,1; 0,0.5; 2,0.5],[ave_R/DD; ave_R/DD; ave_R/DD; ave_R/DD]);
xlabel('$x/D_{C2C}$','FontSize', 18,'Interpreter', 'latex');
ylabel('$y/D_{C2C}$','FontSize', 18,'Interpreter', 'latex');
axis equal;
xlim([0 2]);  ylim([0 1]);




%%

function mu_0 = Get_elastoviscousNum(L, u_0)
% This is to calculate the elasto-viscous number

B = 6.9e-26;  % Bending rigidity
a = 2e-5;  % Diameter of the pillar
mu = 6.1e-3;  % Dynalic viscosity
d = 8e-9; % Diameter of the actin filament

mu_0 = 8 * pi * mu * L^4 * u_0 / (B * a * -log((d/L)^2 * exp(1)));

end


    
    
    
    
