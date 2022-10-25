%%% Calculate the information of rigid fiber cases.

clear; close all; clc;

pathname = uigetdir('F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle', 'Choose a folder');  % input file
Files = dir(fullfile(pathname, 'results with timestamps', '*.mat'));
bg = imread(fullfile(pathname, 'results', "The Pillar.tif"));
the_excel = readcell(fullfile(pathname, 'PoleVaultingCheck.xls'), 'NumHeaderLines', 1);
avi_names = the_excel(:, 1);
im_BW = imbinarize(bg);
LL = regionprops(im_BW, 'Centroid');
Pillar_CoM  = LL.Centroid;
the_pillar_im = imcomplement(imfill(im_BW, 'holes'));
% cmap = cmocean('thermal');
Obj_Mag = 0.63;

im_no = 1:length(Files);
for ii = 1:length(Files)

    load([Files(ii).folder, filesep, Files(ii).name]);
    avi_No = find(cellfun(@(x) any(strcmp(x, append(extractBetween(Files(ii).name, ...
        'trajectory_', '_AABGR'), '.avi'))), avi_names));
    All_data.PoleVaulting(ii) = the_excel{avi_No, 2};
    All_data.Trapping(ii) = the_excel{avi_No, 3};
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
    All_data.timestamps{ii} = Good_case_frm_time - min(Good_case_frm_time);

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

% save([pathname, '_infos_contourL-fromAVG_with_all-Chi_timestamps.mat'], 'All_data');
