%%% Calculate the contact angles and positions

clear; close all; clc;

pathname = uigetdir('F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle', 'Choose a folder');  % input file
Files = dir(fullfile(pathname, 'results with timestamps', '*.mat'));
bg = imread(fullfile(pathname, 'results', "The Pillar (front).tif"));
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
    All_data.acceptability{ii} = the_excel{avi_No, 7};

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
% % %     All_data.initial_y(ii) = ((2048-centroidxy(2, 1)) - Pillar_CoM(2)) * Obj_Mag + 25; % Definition of y0 is the distance between CoM and the edge, so +25.
% % %     All_data.initial_x(ii) = (centroidxy(1, 1) - Pillar_CoM(1)) * Obj_Mag;
% % %     All_data.final_y(ii) = ((2048-centroidxy(2, end)) - Pillar_CoM(2)) * Obj_Mag + 25; % Definition of y0 is the distance between CoM and the edge, so +25.
% % %     All_data.final_x(ii) = (centroidxy(1, end) - Pillar_CoM(1)) * Obj_Mag;
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
    
% % %     % calculate the average speed along x-direction (UNIT: um/s)
% % %     All_data.ave_speed(ii) = ((centroidxy(1, end) - Pillar_CoM(1)) - (centroidxy(1, 1) - Pillar_CoM(1))) * Obj_Mag / (max(Good_case_frm_time) - min(Good_case_frm_time));
% % % 
% % %     % check if the fiber bypass the obtacle tip
% % %     if 2048 - mean(centroidxy(2, and(Pillar_CoM(1)-100 < centroidxy(1, :), centroidxy(1, :) < Pillar_CoM(1)+100))) > Pillar_CoM(2) % if it goes below (in image system).
% % %         bypass_tip = true;
% % %     else
% % %         bypass_tip = false;
% % %     end
% % %     All_data.bypass_tip(ii) = bypass_tip;

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

    for jj = 1:length(Good_case_frm)
        XY = xy.spl{1, Good_case_frm(jj)};
        L_current = sqrt((XY(1,1)-XY(end,1))^2 + (XY(1,2)-XY(end,2))^2);
        if min(pdist2(XY,obs_2d,'euclidean','Smallest',1)) < 8
            Chi_contact = - Chi(jj);

            com_XY = mean(XY, 1); %com_XY(2) = 2048 - com_XY(2);
            L_start = [com_XY(1) - cosd(-Chi_contact) * L_current, com_XY(2) - sind(-Chi_contact)  * L_current];
            L_end = [com_XY(1) + cosd(-Chi_contact) * L_current, com_XY(2) + sind(-Chi_contact)  * L_current];

            P = polyfit(obs_2d(:,1),obs_2d(:,2),1);
            yfit = polyval(P, [900 1200]);

            P_inter = InterX([900 1200; yfit], [L_start; L_end]'); 

            if isempty(P_inter)
                y_contact = nan; 
            else
                y_contact = ((2048-P_inter(2)) - Pillar_CoM(2)) * Obj_Mag + 34; % not so correct because of 'Pillar_CoM' (only front line of the obstacle)
                % Definition of y_c is the distance to the edge, so +25, but 34 is closer to the real value.
            end

%             plot(obs_2d(:,1), obs_2d(:,2)); hold on;
%             plot([L_start(1) L_end(1)], [L_start(2) L_end(2)]); hold on
%             plot(XY(:,1), XY(:,2))
%             close;
            
            break
        end
    end



    try
        All_data.Chi_c(ii) = Chi_contact; All_data.y_c(ii) = y_contact/75;
    catch
        All_data.Chi_c(ii) = nan; All_data.y_c(ii) = nan; 
    end

    clearvars Chi Chi_contact y_contact
end

save([pathname, '_contact_information.mat'], 'All_data');
