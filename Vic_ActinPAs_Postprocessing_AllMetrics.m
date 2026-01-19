%%%
% Calculate bending energy, end-to-end distance, orientation, sphericity,
% curvature, and speed of actin filaments.
%
% Data from Vic_ActinAddPAsInfo.m
%
% Output name format: PlusInfo_trajectory_..._batch1.mat
%%%

clear; close all; clc;

B   = 6.9e-26; % Bending rigidity
mag = 0.1;     % um/pixel

xlsfile = readcell('ForActinPostprocessing.xlsx', 'Sheet', 'Sheet1', 'NumHeaderLines', 1);
NumGroup  = size(xlsfile, 1);
ExpDate   = xlsfile(:, 1);   % Experiment date
storePath = xlsfile(:, 2);   % Path of data to process
Init_U    = xlsfile(:, 8);   % Initial velocity (m/s)
Array_angles = xlsfile(:, 14); % Flow angles

for no_Group = [7 8 13:28]

    the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
    thefiles = dir(fullfile(storePath{no_Group}, '*.mat'));
    Array_angle = Array_angles{no_Group};
    RotMatrix_init = rotz(-Array_angle);
    RotMatrix_init = RotMatrix_init(1:2, 1:2);

    theimgs = dir(['Z:\Experimental Data (EXTRACTED)\Actin Filaments in Porous' ...
        ' Media\', num2str(the_exp_date), '-Actin\AfterAveBGR\*.tif']);

    for file_ind = 1:length(thefiles)

        filename = thefiles(file_ind).name;

        if contains(filename, 'PAsInfoAdded_')

            load(fullfile(thefiles(1).folder, thefiles(file_ind).name));

            pathname       = thefiles(1).folder;
            save_pathname  = strrep(pathname, 'results', 'results_plus_2026');
            if ~isfolder(save_pathname)
                mkdir(save_pathname);
            end
            filename       = thefiles(file_ind).name; 

            % find the corresponding *.tif image
            image_names = struct2cell(theimgs); image_names = image_names(1, :);
            image_ind   = find(cellfun(@(x) contains(x, filename(25:end-11)), image_names));

            spl_Ls = xy.arclen_spl(Good_case_frm);
            L_0    = VicFc_Get_ContourLength(spl_Ls); % filament length

            the_fully_stretched_No = find(abs(spl_Ls - L_0) == min(abs(spl_Ls - L_0)));
            base_fiber_no = Good_case_frm(the_fully_stretched_No);
            base_II = imread(fullfile(theimgs(1).folder, theimgs(image_ind).name), xy(1).frame(base_fiber_no)); % base image

            [intensity_base, ~, ~] = Vic_Fiber_Intensity(xy, the_fully_stretched_No, base_fiber_no, base_II, 30, lzero);

            % Normalize fiber positions w.r.t. pillar lattice
            centers_flip = centers;
            centers_flip(:, 2) = 2048 - centers_flip(:, 2);  % flip to x-y coordinate (NOT image coordinate)
            centers_new = (RotMatrix_init * centers_flip')';  % rotate based on design

            % correct the rotation degree
            sorted_centers_new = sortrows(centers_new, 2);
            sorty = sort(centers_new(:, 2));
            sorty_ind = find(diff(sorty) > 100);
            corrected_angle = zeros(length(sorty_ind)-1, 1);
            for ii = 1:length(sorty_ind)-1
                to_be_fitted = sorted_centers_new(sorty_ind(ii)+1:sorty_ind(ii+1), :);
                fit_linear = fit(to_be_fitted(:, 1), to_be_fitted(:, 2), 'poly1');
                corrected_angle(ii) = atand(fit_linear.p1);
            end
            corrected_angle = mean(corrected_angle);
            RotMatrix_correct = rotz(-Array_angle - corrected_angle);
            RotMatrix_correct = RotMatrix_correct(1:2, 1:2);
            centers_corrected = (RotMatrix_correct * centers_flip')';

            % approximate (grid) pillar center positions
            sortxx = sort(centers_corrected(:, 1));
            sortyy = sort(centers_corrected(:, 2));
            sortxx_ind = find(diff(sortxx) > 100);
            sortyy_ind = find(diff(sortyy) > 100);
            PAs_X = zeros(length(sortxx_ind)+3, 1);
            PAs_Y = zeros(length(sortyy_ind)+3, 1);
            for jj = 1:length(sortxx_ind)-1
                PAs_X(jj+2) = mean(sortxx(sortxx_ind(jj)+1:sortxx_ind(jj+1)));
            end
            for jj = 1:length(sortyy_ind)-1
                PAs_Y(jj+2) = mean(sortyy(sortyy_ind(jj)+1:sortyy_ind(jj+1)));
            end
            % make a grid of pillar CoMs (enlarge the PAs area)
            Ctr2Ctr_x = mean(diff(PAs_X(3:end-2)));
            Ctr2Ctr_y = mean(diff(PAs_Y(3:end-2)));
            PAs_X(2) = PAs_X(3) - Ctr2Ctr_x; PAs_X(1) = PAs_X(2) - Ctr2Ctr_x;
            PAs_X(end-1) = PAs_X(end-2) + Ctr2Ctr_x; PAs_X(end) = PAs_X(end-1) + Ctr2Ctr_x;
            PAs_Y(2) = PAs_Y(3) - Ctr2Ctr_y; PAs_Y(1) = PAs_Y(2) - Ctr2Ctr_y;
            PAs_Y(end-1) = PAs_Y(end-2) + Ctr2Ctr_y; PAs_Y(end) = PAs_Y(end-1) + Ctr2Ctr_y;
            [PAs_X_mesh, PAs_Y_mesh] = meshgrid(PAs_X, PAs_Y);
            PAs_xy = [PAs_X_mesh(:), PAs_Y_mesh(:)];

            nFrm    = size(Good_case_frm, 2);
            CoM_x   = nan(nFrm, 1);
            CoM_y   = nan(nFrm, 1);
            Energy  = nan(nFrm, 1);
            L_ee_norm_belowOne = nan(nFrm, 1);
            L_ee_norm          = nan(nFrm, 1);
            Chi    = nan(nFrm, 1);
            aniso  = nan(nFrm, 1);

            for frm_ind = 1:nFrm

                xy_ind = Good_case_frm(frm_ind); % index of the 'good' case

                crd    = xy.crd{1, xy_ind};
                CoM_xy = xy.centroid{1, xy_ind};

                Gyr = 1/size(crd,1) * [ ...
                    sum((crd(:, 1)-CoM_xy(1)).^2),                sum((crd(:, 1)-CoM_xy(1)) .* (crd(:, 2)-CoM_xy(2))); ...
                    sum((crd(:, 2)-CoM_xy(2)) .* (crd(:, 1)-CoM_xy(1))), sum((crd(:, 2)-CoM_xy(2)).^2) ...
                ];

                [eigenV, eigenD] = eig(Gyr);
                [~, ind] = sort(diag(eigenD));
                Vs = eigenV(:, ind);
                Chi(frm_ind) = atan(Vs(2,2)/Vs(1,2)); % orientation

                Lambda1 = eigenD(2,2); Lambda2 = eigenD(1,1);
                aniso(frm_ind) = 1 - 4*Lambda1*Lambda2/(Lambda1+Lambda2)^2; % sphericity

                II = imread(fullfile(theimgs(1).folder, theimgs(image_ind).name), xy(1).frame(xy_ind)); % load image

                [~, IntensityAll, R2] = Vic_Fiber_Intensity(xy, frm_ind, xy_ind, II, 30, lzero);

                spl_L_current = xy.arclen_spl(xy_ind); % current filament length

                if spl_L_current > 2*L_0
                    CoM_x(frm_ind) = xy.centroid{1, xy_ind}(1);
                    CoM_y(frm_ind) = xy.centroid{1, xy_ind}(2);
                    Energy(frm_ind) = nan;
                else
                    fold_no_max = floor(L_0/spl_L_current) + 1; % max folds based on length

                    intensity_ratio = floor(IntensityAll / intensity_base) - 1;
                    if max(intensity_ratio >= 1)
                        i = max(intensity_ratio);
                        seg_index_tmp = find(intensity_ratio == i);
                        if isscalar(seg_index_tmp)
                            fold_no = i;
                        else
                            fold_no = sum(abs(diff(seg_index_tmp)) > 1) + 1;
                            fold_no = fold_no * i;
                        end
                    else
                        fold_no = 0;
                    end

                    fold_no = min(fold_no, fold_no_max);

                    seglen_spl = xy.seglen_spl{1, xy_ind};
                    CoM_x(frm_ind) = xy.centroid{1, xy_ind}(1);
                    CoM_y(frm_ind) = xy.centroid{1, xy_ind}(2);

                    L_ee_norm(frm_ind) = sqrt((xy.spl{xy_ind}(1,1)-xy.spl{xy_ind}(end,1))^2 ...
                        + (xy.spl{xy_ind}(1,2)-xy.spl{xy_ind}(end,2))^2) / L_0;
                    L_ee_norm_belowOne(frm_ind) = min(L_ee_norm(frm_ind), 1);

                    R2(R2<5) = nan; % remove high-curvature parts (<0.5um radius)

                    Energy(frm_ind) = B / 2 * sum((1./R2(2:end)).^2 .* seglen_spl, 'omitnan') / 1e-7 ...
                        + fold_no * (B / 2 * (1e6)^2 * pi * 1e-6);
                end

                % Curvature for this frame (smoothed)
                spl = xy.spl{1, xy_ind};
                if ~isempty(spl)
                    [~, R2_raw, ~] = curvature(spl);
                    R2_raw(isnan(R2_raw)) = inf;
                    win = max(1, ceil(size(spl,1)/100));
                    RR2 = movmean(R2_raw, win);
                    xy.curvature{1, xy_ind} = 1 ./ RR2; % unit: 1/pixel
                else
                    xy.curvature{1, xy_ind} = [];
                end
            end

            % Normalize fiber CoM based on nearest lower-left pillar
            fiber_CoM_xy = reshape(cell2mat(xy.centroid(1, Good_case_frm)), 2, []);
            Rotated_fiber_CoM_xy = RotMatrix_correct * fiber_CoM_xy;
            [PAs_X_mesh, PAs_Y_mesh] = meshgrid(PAs_X, PAs_Y);
            PAs_xy_complex = complex(PAs_X_mesh(:), PAs_Y_mesh(:));
            Rotated_fiber_CoM_xy_1st = complex(Rotated_fiber_CoM_xy(1, 1), Rotated_fiber_CoM_xy(2, 1));
            real_sign = real(Rotated_fiber_CoM_xy_1st - PAs_xy_complex) > 0;
            imag_sign = imag(Rotated_fiber_CoM_xy_1st - PAs_xy_complex) > 0;
            the_dist_sign = (and(real_sign, imag_sign) - 0.5) * 2;
            dist = the_dist_sign .* abs(Rotated_fiber_CoM_xy_1st - PAs_xy_complex);
            positive_indices = find(dist > 0);
            min_positive_value = min(dist(positive_indices));
            min_positive_position = find(dist == min_positive_value);
            ref_PAs_xy = PAs_xy_complex(min_positive_position);
            Rotated_fiber_CoM_xy_normalized = [Rotated_fiber_CoM_xy(1, :) - real(ref_PAs_xy); ...
                Rotated_fiber_CoM_xy(2, :) - imag(ref_PAs_xy)];

            % Normalize fiber coordinates in lattice
            for frm_ind = 1:size(Good_case_frm, 2)
                xy_ind = Good_case_frm(frm_ind);
                spl = xy.spl{1, xy_ind};
                Rotated_spl = (RotMatrix_correct * spl')';
                D = pdist2(PAs_xy, Rotated_spl);
                PAs_x_repeat = repmat(PAs_X_mesh(:), 1, size(Rotated_spl, 1));
                PAs_y_repeat = repmat(PAs_Y_mesh(:), 1, size(Rotated_spl, 1));

                spl_x_repeat = repmat(Rotated_spl(:, 1)', size(PAs_xy, 1), 1);
                spl_y_repeat = repmat(Rotated_spl(:, 2)', size(PAs_xy, 1), 1);
                dist_x = spl_x_repeat - PAs_x_repeat;
                dist_y = spl_y_repeat - PAs_y_repeat;
                real_sign_spl = dist_x > 0;
                imag_sign_spl = dist_y > 0;
                the_dist_sign_spl = (and(real_sign_spl, imag_sign_spl) - 0.5) * 2;
                dist_signed = the_dist_sign_spl .* D;
                dist_signed(dist_signed < 0) = nan;
                [ref_x_ind, ref_y_ind] = find(dist_signed == min(dist_signed));

                fiber_xy_in_lattice = [diag(dist_x(ref_x_ind, ref_y_ind)) / Ctr2Ctr_x, ...
                    diag(dist_y(ref_x_ind, ref_y_ind)) / Ctr2Ctr_y];
                xy.fiber_xy_in_lattice{1, xy_ind} = fiber_xy_in_lattice;
            end

            % Normalize fiber CoM in lattice
            Rotated_CoM_xy = (RotMatrix_correct * fiber_CoM_xy)';
            D_CoM = pdist2(PAs_xy, Rotated_CoM_xy);
            PAs_x_repeat_CoM = repmat(PAs_X_mesh(:), 1, size(Rotated_CoM_xy, 1));
            PAs_y_repeat_CoM = repmat(PAs_Y_mesh(:), 1, size(Rotated_CoM_xy, 1));
            CoM_x_repeat = repmat(Rotated_CoM_xy(:, 1)', size(PAs_xy, 1), 1);
            CoM_y_repeat = repmat(Rotated_CoM_xy(:, 2)', size(PAs_xy, 1), 1);
            dist_CoM_x = CoM_x_repeat - PAs_x_repeat_CoM;
            dist_CoM_y = CoM_y_repeat - PAs_y_repeat_CoM;
            real_sign_CoM = dist_CoM_x > 0;
            imag_sign_CoM = dist_CoM_y > 0;
            the_dist_sign_CoM = (and(real_sign_CoM, imag_sign_CoM) - 0.5) * 2;
            dist_CoM = the_dist_sign_CoM .* D_CoM;
            dist_CoM(dist_CoM < 0) = nan;
            [ref_x_ind_CoM, ref_y_ind_CoM] = find(dist_CoM == min(dist_CoM));
            fiber_CoM_xy_in_lattice = [diag(dist_CoM_x(ref_x_ind_CoM, ref_y_ind_CoM)) / Ctr2Ctr_x, ...
                diag(dist_CoM_y(ref_x_ind_CoM, ref_y_ind_CoM)) / Ctr2Ctr_y];

            % Instantaneous velocity (normalized by Init_U)
            if contains(filename, 'interval4')
                dt = 0.08;
            elseif contains(filename, 'interval2')
                dt = 0.04;
            else
                dt = 0.02;
            end
            Instan_V = (sqrt(diff(CoM_x).^2 + diff(CoM_y).^2))' * mag * 1e-6 ./ ...
                (diff(xy(1).frame(Good_case_frm)) * dt) / Init_U{no_Group};

            % Build save list including normalization outputs
            vars_to_save = {'centers','circleMask','CoM_x','CoM_y','Energy','Good_case_frm', ...
                'L_ee_norm','L_ee_norm_belowOne','lzero','metric','radii','xy','Chi','aniso', ...
                'Instan_V','Rotated_fiber_CoM_xy_normalized','RotMatrix_correct','Ctr2Ctr_x', ...
                'Ctr2Ctr_y','fiber_CoM_xy_in_lattice'};
            if exist('fiber_xy_in_lattice','var')
                vars_to_save{end+1} = 'fiber_xy_in_lattice'; 
            end

            save([save_pathname, filesep, 'PlusInfo_', filename(14: end-4), '.mat'], vars_to_save{:})

            clearvars CoM_x CoM_y Energy xy centers radii L_ee_norm_belowOne L_ee_norm Chi aniso Instan_V ...
                      Rotated_fiber_CoM_xy_normalized RotMatrix_correct Ctr2Ctr_x Ctr2Ctr_y fiber_CoM_xy_in_lattice fiber_xy_in_lattice
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (Vic_Get_ave_cutExtreme and Vic_Fiber_Intensity unchanged)
function L_0 = Vic_Get_ave_cutExtreme(spl_Ls, threshold)
    L_0_coarse = mean(spl_Ls, 'omitnan');
    spl_Ls(spl_Ls < L_0_coarse*(1-threshold)) = [];
    spl_Ls(spl_Ls > L_0_coarse*(1+threshold)) = [];
    if isempty(spl_Ls); L_0 = L_0_coarse; return; end
    L_0 = mean(spl_Ls, 'omitnan');
    while L_0_coarse ~= L_0 && ~isnan(L_0)
        L_0_coarse = L_0;
        spl_Ls(spl_Ls < L_0_coarse*(1-threshold)) = [];
        spl_Ls(spl_Ls > L_0_coarse*(1+threshold)) = [];
        if isempty(spl_Ls); L_0 = L_0_coarse; return; end
        L_0 = mean(spl_Ls, 'omitnan');
    end
end

function [intensity_0, IntensityAll, R] = Vic_Fiber_Intensity(xy, frm_ind, xy_ind, II, cut_length, lzero)
    spl = xy(1).spl{xy_ind};
    spl(:, 1) = spl(:, 1) + lzero(frm_ind);
    spl(:, 2) = size(II, 2) - spl(:, 2) - lzero(frm_ind);
    [~,R,K] = curvature(spl);
    IntensityAll = zeros(size(spl,1), 1);
    for seg_ind = 1:size(spl,1)
        P_0 = spl(seg_ind, :); vec_0 = K(seg_ind, :);
        if ~isnan(vec_0(1))
            P_start = P_0 + cut_length/2 * vec_0/norm(vec_0);
            P_end   = P_0 - cut_length/2 * vec_0/norm(vec_0);
            The_cut = round([linspace(P_start(1), P_end(1), cut_length)', ...
                             linspace(P_start(2), P_end(2), cut_length)']);
        else
            try
                vec_1 = [spl(seg_ind, 1) - spl(seg_ind+1, 1), ...
                         spl(seg_ind, 2) - spl(seg_ind+1, 2)];
            catch
                vec_1 = [spl(seg_ind, 1) - spl(seg_ind-1, 1), ...
                         spl(seg_ind, 2) - spl(seg_ind-1, 2)];
            end
            if isempty(find(vec_1 == 0, 1))
                vec_0 = [1, -vec_1(1)/vec_1(2)]; vec_0 = vec_0/norm(vec_0);
            else
                vec_0 = zeros(1, 2); vec_0(vec_1 == 0) = 1;
            end
            P_start = P_0 + cut_length/2 * vec_0;
            P_end   = P_0 - cut_length/2 * vec_0;
            The_cut = round([linspace(P_start(1), P_end(1), cut_length)', ...
                             linspace(P_start(2), P_end(2), cut_length)']);
        end

        try
            cut_intensity = II(sub2ind([size(II, 1), size(II, 2)], The_cut(:,2), The_cut(:,1)));
            XX = (1:cut_length)';
            try
                intensity_fit = fit(XX,double(cut_intensity),'gauss1');
                intensity_fit_fun = @(x) intensity_fit.a1.*exp(-((x-intensity_fit.b1)./intensity_fit.c1).^2);
                if intensity_fit.b1 > 1 && intensity_fit.b1 < cut_length
                    IntensityAll(seg_ind) = integral(intensity_fit_fun,-Inf,Inf);
                else
                    IntensityAll(seg_ind) = nan;
                end
            catch
                IntensityAll(seg_ind) = nan;
            end
        catch
            IntensityAll(seg_ind) = nan;
        end
    end
    intensity_0 = Vic_Get_ave_cutExtreme(IntensityAll, 0.2);
end
