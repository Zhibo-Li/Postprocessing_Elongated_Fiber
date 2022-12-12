%%% Calculate the information of rigid fiber cases.

clear; close all; clc;

pathname = uigetdir('F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle', 'Choose a folder');  % input file
Files = dir(fullfile(pathname, 'results with timestamps', '*.mat'));
bg = imread(fullfile(pathname, 'results', "The Pillar.tif"));
the_excel = readcell(fullfile(pathname, 'PoleVaultingCheck.xls'), 'NumHeaderLines', 1);
avi_names = the_excel(:, 1);
im_BW = imbinarize(bg);
[row, col] = find(im_BW); obs_2d = [col, 2048-row];
LL = regionprops(im_BW, 'Centroid');
Pillar_CoM  = LL.Centroid;
the_pillar_im = imcomplement(imfill(im_BW, 'holes'));
% cmap = cmocean('thermal');
Obj_Mag = 0.63;

% load LBM simulation data 
load('D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\LBM data\LBM_midplane.mat');

im_no = 1:length(Files);
for ii = 1:length(Files)

    load([Files(ii).folder, filesep, Files(ii).name]);
    avi_No = find(cellfun(@(x) any(strcmp(x, append(extractBetween(Files(ii).name, ...
        'trajectory_', '_AABGR'), '.avi'))), avi_names));
    All_data.PoleVaulting(ii) = the_excel{avi_No, 2};
    ifTrapping = the_excel{avi_No, 3};
    All_data.Trapping(ii) = ifTrapping;
    All_data.Sliding(ii) = the_excel{avi_No, 4};
    All_data.ApexVaulting(ii) = the_excel{avi_No, 5};

%     %%% Plot trajectories
%     the_time = Good_case_frm_time * 100; % Here * 100 for convenience.
%     the_time = round(the_time - min(the_time));
%     
%     figure('color', 'w'); set(gcf, 'Position', [300 300 600 100]);
%     imshow(the_pillar_im(801:1200, :)); hold on
%     %%% define the range of loop indicator k based on the colormap
% %     if max(the_time) > 255
% %         k_range = find(the_time > 255, 1) - 1;
% %         for k =1:min(k_range, length(Good_case_frm))
% %             plot(xy.spl{Good_case_frm(k)}(:,1),2048-xy.spl{Good_case_frm(k)}(:,2)-800,'Color',cmap(the_time(k)+1,:),'LineWidth',1.5)
% %             hold on
% %         end
% %     else
% %         for k =1:length(Good_case_frm)
% %             plot(xy.spl{Good_case_frm(k)}(:,1),2048-xy.spl{Good_case_frm(k)}(:,2)-800,'Color',cmap(the_time(k)+1,:),'LineWidth',1.5)
% %             hold on
% %         end
% %     end
%     %%% draw trajectories with looping the colormap
%     for k =1:length(Good_case_frm)
%         plot(xy.spl{Good_case_frm(k)}(:,1),2048-xy.spl{Good_case_frm(k)}(:,2)-800,'Color',cmap(mod(the_time(k)+1, 256)+1,:),'LineWidth',1.5)
%         hold on
%     end
%     axis equal
%     hold off
% %     hcb=colorbar;
% %     title(hcb,'$Time\ (s)$','FontSize', 16,'Interpreter', 'latex');
%     f=gcf;
%     exportgraphics(f,[pathname,'\results with timestamps\',Files(ii).name,'_looping_color.png'],'Resolution',100)
%     close

    % add the names
    All_data.filename{ii} = Files(ii).name;    

    % add the timestamps
    timestamps = Good_case_frm_time - min(Good_case_frm_time);
    All_data.timestamps{ii} = timestamps;

    % add the contour lengths 
    sorted_lengths = sort(xy.arclen_spl(Good_case_frm));
    All_data.contour_length(ii) = mean(sorted_lengths) * Obj_Mag;  % Average to get the contour length (UNIT: um).

    % add the position information
    centroidxy = reshape(cell2mat(xy.centroid),2,numel(xy.centroid));
    centroidxy = centroidxy(:, Good_case_frm);
    All_data.delta_y(ii) = (centroidxy(2, 1) - centroidxy(2,end)) * Obj_Mag;
    All_data.initial_y(ii) = ((2048-centroidxy(2, 1)) - Pillar_CoM(2)) * Obj_Mag + 25; % Definition of y0 is the distance between CoM and the edge, so +25.
    All_data.initial_x(ii) = (centroidxy(1, 1) - Pillar_CoM(1)) * Obj_Mag;
    All_data.final_y(ii) = ((2048-centroidxy(2, end)) - Pillar_CoM(2)) * Obj_Mag + 25; % Definition of y0 is the distance between CoM and the edge, so +25.
    All_data.final_x(ii) = (centroidxy(1, end) - Pillar_CoM(1)) * Obj_Mag;
    All_data.CoM{ii} = centroidxy;
    %     hold on; plot(centroidxy(1,:), 2048-centroidxy(2,:)-800, "Color",'k', 'LineStyle','--');

    % calculate the interaction index
    dt = diff(timestamps); % [time(i+1) - time(i)]
    dx = diff(centroidxy(1, :))'; % [x(i+1) - x(i)]
    dx_dt = movmean(dx, 7) ./ movmean(dt, 7) * Obj_Mag; % [chi(i+1) - chi(i)] / [time(i+1) - time(i)]. 
    num_up_obs = sum(centroidxy(1, :) < 300); % set range (without the perturbation of the obstacle) for average speed calculation (in pixel).
    speed_upstream = mean(dx_dt(1:num_up_obs)); % average speed (without the perturbation of the obstacle)
    num_down_obs = sum(centroidxy(1, :) > 1749); % set range (without the perturbation of the obstacle) for average speed calculation (in pixel).
    if ifTrapping == 1
        speed_downstream = speed_upstream;
    else
        speed_downstream = mean(dx_dt(end-num_down_obs+1:end)); % average speed (without the perturbation of the obstacle)
    end
    U0 = 0.5 * (speed_upstream + speed_downstream); % U0: average of speed_upstream and speed_downstream
    speed_ave_all = mean(dx_dt); % average all instantaneous speeds
    All_data.speed_ave_all(ii) = speed_ave_all;
    All_data.speed_upstream(ii) = speed_upstream;
    All_data.speed_downstream(ii) = speed_downstream;
    % interaction index 1 (defines as 1-speed_ave_all/U0)
    interaction1 = 1-speed_ave_all/U0;
    All_data.interaction1(ii) = interaction1;

    % interaction index 2 (defines as max(abs(U0-U(t)))/U0)
    interaction2 = max(abs(dx_dt-U0))/U0; 
    All_data.interaction2(ii) = interaction2;
    
    % interaction index 3 (defines as contact probability -- varying layer thickness)
    circle_intersec = 0; line_intersec = 0;
    for jj = 1:size(centroidxy, 2)

        XY = xy.spl{1, Good_case_frm(jj)}; 
        Half_L = sqrt((XY(1,1)-XY(end,1))^2 + (XY(1,2)-XY(end,2))^2) / 2 + 8;

        theta = linspace(0,2*pi,300);
        circle_x = Half_L*cos(theta); circle_x = circle_x + centroidxy(1, jj);
        circle_y = Half_L*sin(theta); circle_y = circle_y + centroidxy(2, jj);
        [intersec_cir_x, intersec_cir_y] = polyxpoly(obs_2d(:,1), obs_2d(:,2), circle_x, circle_y);
        if ~isempty(intersec_cir_y)
            circle_intersec = circle_intersec + 1;
            if min(pdist2(XY,obs_2d,'euclidean','Smallest',1)) < 8
                line_intersec = line_intersec + 1;
            end
        end
    end
    if circle_intersec == 0
        interaction3 = 0;
    else
        interaction3 = line_intersec / circle_intersec;
    end
    All_data.interaction3(ii) = interaction3;

    % interaction index 4 (defines as the ratio of 'integral of the speed-U0 difference' and 'integral of the U0')
    I_a = sum(U0.*diff(timestamps(1:end-1))); % integral of the U0.
    I_b = sum(abs(dx_dt(1:end-1)-U0).*diff(timestamps(1:end-1))); % integral of the speed-U0 difference.
    All_data.interaction4(ii) = I_b/I_a;

    % interaction index 5 (defines as the normalized integral of speed difference compared to the streamline)
    relative_x = (centroidxy(1, :) - Pillar_CoM(1)) * Obj_Mag; 
    relative_y = (2048 - centroidxy(2, :) - Pillar_CoM(1)) * Obj_Mag; % relative x & y to the obstacle center.
    obs_simu_ctr = [XX(1, 151), YY(51, 1)]; % The obstacle center in the simulation.
    Ux_simu = zeros(size(relative_x, 2), 1); Uy_simu = Ux_simu;
    for kk = 1:size(relative_x, 2)
        XY_simu = [relative_x(kk) + obs_simu_ctr(1); relative_y(kk) + obs_simu_ctr(2)];

        Ux_simu(kk) = interp2(XX, YY, UX_midplane, XY_simu(1), XY_simu(2));
        Uy_simu(kk) = interp2(XX, YY, UY_midplane, XY_simu(1), XY_simu(2));
    end
    Ux_simu = Ux_simu / Ux_simu(1) * U0; Uy_simu = Uy_simu / Uy_simu(1) * U0; % scale the simulation velocity.
    All_data.interaction5(ii) = sum(abs(Ux_simu(1:end-1) - dx_dt) .* diff(timestamps)) / (timestamps(end) * U0);

    % interaction index 6 (defines as contact probability -- fixed layer thickness of 40um)
    circle_intersec = 0; line_intersec = 0;
    for jj = 1:size(centroidxy, 2)

        XY = xy.spl{1, Good_case_frm(jj)};
        
        if min(pdist2(centroidxy(:, jj)',obs_2d,'euclidean','Smallest',1)) < (40 / Obj_Mag)
            circle_intersec = circle_intersec + 1;
            if min(pdist2(XY,obs_2d,'euclidean','Smallest',1)) < 8
                line_intersec = line_intersec + 1;
            end
        end
    end
    if circle_intersec == 0
        interaction6 = 0;
    else
        interaction6 = line_intersec / circle_intersec;
    end
    All_data.interaction6(ii) = interaction6;

    % calculate the average speed along x-direction (UNIT: um/s)
    All_data.ave_speed(ii) = ((centroidxy(1, end) - Pillar_CoM(1)) - (centroidxy(1, 1) - Pillar_CoM(1))) * Obj_Mag / (max(Good_case_frm_time) - min(Good_case_frm_time));

    % check if the fiber bypass the obtacle tip
    if 2048 - mean(centroidxy(2, and(Pillar_CoM(1)-100 < centroidxy(1, :), centroidxy(1, :) < Pillar_CoM(1)+100))) > Pillar_CoM(2) % if it goes below (in image system).
        bypass_tip = true;
    else
        bypass_tip = false;
    end
    All_data.bypass_tip(ii) = bypass_tip;

    % calculate the initial angle Chi_0
    crd = xy.crd{1,Good_case_frm(1)};
    centroidxy_k = xy.centroid{1,Good_case_frm(1)};
    Gyr = 1/size(crd,1) * [sum((crd(:, 1)-centroidxy_k(1)).^2),  sum((crd(:, 1)-centroidxy_k(1)) .* (crd(:, 2)-centroidxy_k(2)));
        sum((crd(:, 2)-centroidxy_k(2)) .* (crd(:, 1)-centroidxy_k(1))), sum((crd(:, 2)-centroidxy_k(2)).^2)];

    [eigenV,eigenD] = eig(Gyr);
    [d,ind] = sort(diag(eigenD));
    Ds = eigenD(ind,ind);
    Vs = eigenV(:,ind);
    All_data.Chi_0(ii)  = atan(Vs(2,2)/Vs(1,2))/pi*180;

    % calculate all the angles
    for jj = 1:length(Good_case_frm)
        crd = xy.crd{1,Good_case_frm(jj)};
        centroidxy_k = xy.centroid{1,Good_case_frm(jj)};
        Gyr = 1/size(crd,1) * [sum((crd(:, 1)-centroidxy_k(1)).^2),  sum((crd(:, 1)-centroidxy_k(1)) .* (crd(:, 2)-centroidxy_k(2)));
            sum((crd(:, 2)-centroidxy_k(2)) .* (crd(:, 1)-centroidxy_k(1))), sum((crd(:, 2)-centroidxy_k(2)).^2)];

        [eigenV,eigenD] = eig(Gyr);
        [d,ind] = sort(diag(eigenD));
        Ds = eigenD(ind,ind);
        Vs = eigenV(:,ind);
        Chi(jj) = atan(Vs(2,2)/Vs(1,2))/pi*180;
    end
    All_data.Chi{ii}  = Chi;
    clearvars Chi

end

save([pathname, '_information_full.mat'], 'All_data');
