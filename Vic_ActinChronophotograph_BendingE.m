clear; close all; clc;

% % % % The .mat file where stores your results.
% [filename, pathname]=uigetfile({['F:\Processing & Results\Actin Filaments in Porous Media\' ...
%     '20210914-Actin\results\Group_1\*.mat']}, 'Choose a *.mat to be processed');  % input file
% load([pathname, filename])

% % % The .mat file where stores the PAs information.
load(['Z:\Processing & Results\Actin Filaments in Porous Media\20220216-Actin' ...
    '\circlesforPAs20220216_T10.mat'])
centers(:, 2) = 2048 - centers(:, 2);

% % % The .mat files where stores your results.
thefiles = dir(['Z:\Processing & Results\Actin Filaments in Porous Media' ...
    '\20220216-Actin\results\*.mat']);

% % % The .tif files where stores your results.
theimgs = dir(['Z:\Experimental Data (EXTRACTED)\Actin Filaments in Porous Media' ...
    '\20220216-Actin\AfterAveBGR\*.tif']);

% % % The .mat file where stores the PAs information.
% % % load(['D:\Dropbox\PROCESS remotely\20230118\New calculations\20210914-Actin' ...
% % %     '\results\circlesforPAs20210914_G10.mat'])
% % % centers(:, 2) = 2048 - centers(:, 2);
% % % 
% % % % % % The .mat files where stores your results.
% % % thefiles = dir(['D:\Dropbox\PROCESS remotely\20230118\New calculations' ...
% % %     '\20210914-Actin\results\Group_1\*.mat']);
% % % 
% % % % % % The .tif files where stores your results.
% % % theimgs = dir('D:\Dropbox\PROCESS remotely\20230118\*.tif');

B = 6.9e-26;  % Bending rigidity
% base_II = imread(['F:\Experimental Data (EXTRACTED)\Actin Filaments in Porous' ...
%     ' Media\20220216-Actin\The masks\T10_23.08.21_upperone_1.tif']);
% II = base_II;

for file_ind = 27:49%1:length(thefiles)

    load(fullfile(thefiles(1).folder, thefiles(file_ind).name));

    pathname = thefiles(1).folder;
    filename = thefiles(file_ind).name

    % find the cooresponding *.tif image
    image_names = struct2cell(theimgs); image_names = image_names(1, :);
    image_ind = find(cellfun(@(x) contains(x, filename(12:end-11)), image_names));

    figure('color', 'w'); set(gcf, 'Position', [100 300 1600 300]);
    cmap = colormap('jet');

%     Good_case_frm = Good_case;
    Good_case_frm = find(ismember(xy(1).frame, Good_case));
    spl_Ls = xy.arclen_spl(Good_case_frm);
    L_0 = Vic_Get_ave_cutExtreme(spl_Ls, 0.2); % get the filament length
%     L_0 = Vic_Get_L_histo(spl_Ls); % get the filament length (not good)

    the_fully_stretched_No = find(abs(spl_Ls - L_0) == min(abs(spl_Ls - L_0))); 
    base_fiber_no = Good_case_frm(the_fully_stretched_No);
    base_II = imread(fullfile(theimgs(1).folder, theimgs(image_ind).name), xy(1).frame(base_fiber_no)); % load the image
    lzero = 0; % careful (there is another line.)
%     try
%         lzero = max(lobject,ceil(5*lnoise)); % important!
%     catch
%         lzero = max(prmt(base_fiber_no).lobject,ceil(5*prmt(base_fiber_no).lnoise));
%     end
    [intensity_base, ~, ~] = Vic_Fiber_Intensity(xy, base_fiber_no, base_II, 30, lzero);

    CoM_x = zeros(size(Good_case_frm,2), 1); Energy = zeros(size(Good_case_frm,2), 1);
    for frm_ind = 1:size(Good_case_frm,2)

        xy_ind = Good_case_frm(frm_ind);% index of the 'good' cases
        II = imread(fullfile(theimgs(1).folder, theimgs(image_ind).name), xy(1).frame(xy_ind)); % load the image
        lzero = 0; % careful (there is another line.)
%         try
%             lzero = max(lobject,ceil(5*lnoise)); % important!
%         catch
%             lzero = max(prmt(xy_ind).lobject,ceil(5*prmt(xy_ind).lnoise));
%         end
        %%%%%%%%%%%%%% plot fiber snapshots %%%%%%%%%%%%%%%%%%%%%%
        if frm_ind == 1
            plot(xy.spl{xy_ind}(:,1)+lzero,xy.spl{xy_ind}(:,2)+lzero, 'LineWidth', 2);
            addaxislabel(1,'y (pixel)');
        elseif mod(frm_ind, 3) == 0
            addaxisplot(xy.spl{xy_ind}(:,1)+lzero,xy.spl{xy_ind}(:,2)+lzero, 1, 'color', cmap(mod(frm_ind*32, 255)+1, :), 'LineWidth', 2);
        end
        hold on
        %%%%%%%%%%%%%% plot fiber snapshots %%%%%%%%%%%%%%%%%%%%%%

        [~, IntensityAll,R2] = Vic_Fiber_Intensity(xy, xy_ind, II, 30, lzero);

        spl_L_current = xy.arclen_spl(xy_ind); % fiber length in the current processing frame

        if spl_L_current > 2*L_0 % discard the too-long fiber
            CoM_x(frm_ind) = xy.centroid{1, xy_ind}(1); % x-position of fiber center-of-mass
            Energy(frm_ind) = nan;
        else
            fold_no_max = floor(L_0/spl_L_current) + 1; % maximum folder number based on the fiber length

            intensity_ratio = floor(IntensityAll / intensity_base) - 1;
            if max(intensity_ratio >= 1)
                i = max(intensity_ratio);
                seg_index_tmp = find(intensity_ratio == i);
                if length(seg_index_tmp) == 1
                    fold_no = i;
                else
                    fold_no = sum(abs(diff(seg_index_tmp)) > 1) + 1;
                    fold_no = fold_no * i; % the number of folders on the filament
                end

            else
                fold_no = 0;
            end

            fold_no = min(fold_no, fold_no_max);

            seglen_spl = xy.seglen_spl{1, xy_ind}; % segment lengths
            CoM_x(frm_ind) = xy.centroid{1, xy_ind}(1); % x-position of fiber center-of-mass

            R2(R2<5) = nan; % remove the too-high curvature part (here, the allowed minimum curvature radius is 0.5um).

            Energy(frm_ind) = B / 2 * sum((1./R2(2:end)).^2 .* seglen_spl, ...
                'omitnan') / 1e-7 + fold_no * (B / 2 * (1e6)^2 * pi * 1e-6);  
            % 1e-7 is the scale m/pixel;
        end
    end
    xlim([0 2050])
    xlabel('x (pixel)')
    axis equal; hold on

    addaxis(CoM_x+lzero, Energy, '*r', 'LineStyle','none', 'MarkerSize', 7);
    addaxislabel(2, 'Bending energy (J)');
    viscircles(centers, radii,'LineStyle','--', 'LineWidth', 0.5, 'Color', 'k'); hold on
    
    f=gcf;
    exportgraphics(f,[pathname, filesep, filename(1: end-4), '_BendingEnergy_new.png'],'Resolution',100)

    close all
    clearvars CoM_x Energy xy
end






function L_0 = Vic_Get_ave_cutExtreme(spl_Ls, threshold)
% the function is used to calculte the contour lenght of the filament based
% on the average without the extrame values.
% 
% threshold -- threshold to truncate the extreme values 

    L_0_coarse = mean(spl_Ls, 'omitnan'); % roughly estimated filament length
    spl_Ls(spl_Ls < L_0_coarse*(1-threshold)) = []; 
    spl_Ls(spl_Ls > L_0_coarse*(1+threshold)) = []; % remove the extreme value
    if isempty(spl_Ls)
        L_0 = L_0_coarse; 
        return 
    end
    L_0 = mean(spl_Ls, 'omitnan');

    while L_0_coarse ~= L_0 && ~isnan(L_0)
        L_0_coarse = L_0;
        spl_Ls(spl_Ls < L_0_coarse*(1-threshold)) = [];
        spl_Ls(spl_Ls > L_0_coarse*(1+threshold)) = []; % remove the extreme value
        if isempty(spl_Ls)
            L_0 = L_0_coarse;
            return
        end
        L_0 = mean(spl_Ls, 'omitnan');
    end
end

function L_0 = Vic_Get_L_histo(spl_Ls)
% the function is used to calculte the contour lenght of the filament based
% on histogram.

    [N,edges] = histcounts(spl_Ls, 5); % here choose 10 bins
    L_0 = mean(edges(find(N == max(N)) + 1) + edges(N == max(N)));
end


function [intensity_0, IntensityAll, R] = Vic_Fiber_Intensity(xy, xy_ind, II, cut_length, lzero)
% get the intensity along the fiber
% xy -- the struct contains the reconstructed information
% xy_ind -- current index of xy
% II -- the corresponding raw image
% cut_length -- the length of normal cut (should be even number)
% lzero -- a useless parameter (but necessary).
%
% IntensityAll -- the integral along the centerline of the filament
% intensity_0 -- the average or common value along the centerline of the filament
% R2 -- Radius of curvature along the fiber

spl = xy(1).spl{xy_ind}; % the x-y coordinates of the spline (in 'plot' coordinate)
spl(:, 1) = spl(:, 1) + lzero; % lzero is very important !!!
spl(:, 2) = size(II, 2) - spl(:, 2) - lzero; % convert the x-y coordinates to 'image' coordinate

[~,R,K] = curvature(spl); % calculate the curvature
% RR2 = movmean(R2,ceil(size(spl,1)/100)); % ???????????? needed???????

IntensityAll = zeros(size(spl,1), 1);
for seg_ind = 1:size(spl,1)

    P_0 = spl(seg_ind, :); vec_0 = K(seg_ind, :); % given point and normal vector
    if ~isnan(vec_0(1))
        P_start = P_0 + cut_length/2 * vec_0/norm(vec_0);
        P_end = P_0 - cut_length/2 * vec_0/norm(vec_0); % normal to the spline
        The_cut = round([linspace(P_start(1), P_end(1), cut_length)', ...
            linspace(P_start(2), P_end(2), cut_length)']); % the normal cut
    else
        try
            vec_1 = [spl(seg_ind, 1) - spl(seg_ind+1, 1), ...
                spl(seg_ind, 2) - spl(seg_ind+1, 2)]; % tangent vector
        catch
            vec_1 = [spl(seg_ind, 1) - spl(seg_ind-1, 1), ...
                spl(seg_ind, 2) - spl(seg_ind-1, 2)]; % tangent vector (if the end-point on the fiber)
        end
        if isempty(find(vec_1 == 0, 1))
            vec_0 = [1, -vec_1(1)/vec_1(2)]; vec_0 = vec_0/norm(vec_0);
        else
            vec_0 = zeros(1, 2); vec_0(vec_1 == 0) = 1; % normal vector
        end
        P_start = P_0 + cut_length/2 * vec_0;
        P_end = P_0 - cut_length/2 * vec_0; % normal to the spline
        The_cut = round([linspace(P_start(1), P_end(1), cut_length)', ...
            linspace(P_start(2), P_end(2), cut_length)']); % the normal cut
    end

%     imshow(II, []);hold on
%     plot(xy.spl{xy_ind}(:,1)+lzero,size(II, 2)-xy.spl{xy_ind}(:,2)-lzero, 'LineWidth', 2);hold on
%     plot(The_cut(:,1),The_cut(:,2), 'LineWidth', 2);hold off

    try
        cut_intensity = II(sub2ind([size(II, 1), size(II, 2)], ...
            The_cut(:,2), The_cut(:,1))); % the intensity along the normal cut

        XX = (1:cut_length)'; % unit: pixel
        try
            intensity_fit = fit(XX,double(cut_intensity),'gauss1'); % sometimes it fails
%             figure; plot(intensity_fit, XX, cut_intensity)
            intensity_fit_fun = @(x) intensity_fit.a1.*exp(-((x-intensity_fit.b1)./intensity_fit.c1).^2);
            if intensity_fit.b1 > 1 && intensity_fit.b1 < cut_length
                IntensityAll(seg_ind) = integral(intensity_fit_fun,-Inf,Inf);
%                 the integral of the intensity along the cut line
            else
                IntensityAll(seg_ind) = nan;
%                 remove the unreasonable integral value
            end
        catch
            IntensityAll(seg_ind) = nan;
        end
    catch
        IntensityAll(seg_ind) = nan; % set to be nan if The_cut is too long and out of the image
    end

end

% figure; plot(IntensityAll)
intensity_0 = Vic_Get_ave_cutExtreme(IntensityAll, 0.2);

end