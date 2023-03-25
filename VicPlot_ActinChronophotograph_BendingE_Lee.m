%%%% Plot the chronophotograoh, bending energy and Lee of actin filaments.
%%%% data from VicPlot_ActinChronophotograph_BendingE_Lee.m (the second section)
%%%% data name format: trajectory_..._batch1_Info.mat

clear; close all; clc;

set(0,'DefaultAxesFontSize',14);
set(0,'defaulttextfontsize',14);
% set(0,'defaultAxesFontName', 'times new roman');
% set(0,'defaultTextFontName', 'times new roman');
set(0,'defaultAxesFontName', 'Arial');
set(0,'defaultTextFontName', 'Arial');
set(0,'defaulttextInterpreter','latex') 

mag = 0.1; % um/pixel

xlsfile = readcell('ForActinPostprocessing.xlsx','Sheet','Sheet1','NumHeaderLines',1);
% This is the file that contains all the information about the later processing (in sheet 1).

NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.

for no_Group = [7 8 13:25]

    the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
    thefiles = dir(fullfile(storePath{no_Group},'*.mat'));

    for file_ind = 1:length(thefiles)

        filename = thefiles(file_ind).name;

        if contains(filename, '_Info') 

            load(fullfile(thefiles(1).folder, thefiles(file_ind).name));

            pathname = thefiles(1).folder;
            filename = thefiles(file_ind).name

            figure('color', 'w'); set(gcf, 'Position', [100 300 1600 400]);
            cmap = colormap('jet');

            for frm_ind = 1:size(Good_case_frm,2)

                xy_ind = Good_case_frm(frm_ind);% index of the 'good' cases

                %%%%%%%%%%%%%% plot fiber snapshots %%%%%%%%%%%%%%%%%%%%%%
                if frm_ind == 1
                    plot((xy.spl{xy_ind}(:,1)+lzero(frm_ind))*mag, ...
                        (xy.spl{xy_ind}(:,2)+lzero(frm_ind))*mag, 'LineWidth', 2);
%                     addaxislabel(1,'y (pixel)');
                elseif mod(frm_ind, 3) == 0
                    addaxisplot((xy.spl{xy_ind}(:,1)+lzero(frm_ind))*mag, ...
                        (xy.spl{xy_ind}(:,2)+lzero(frm_ind))*mag, 1, 'color', cmap(mod(frm_ind*32, 255)+1, :), 'LineWidth', 2);
                end
                hold on

            end
            set(gca,'box','off'); set(gca,'ytick',[]);
            xlim([0 2050]*mag); 
            xlabel('$x\ (\mu{m})$')
            axis equal; hold on

%             addaxis((CoM_x+lzero)*mag, L_ee_norm, [0 1],'*r', 'LineStyle','none', 'MarkerSize', 10);
            addaxis((CoM_x+lzero)*mag, L_ee_norm,'*r', 'LineStyle','none', 'MarkerSize', 10);
            addaxislabel(2, '$L_{ee}/L_0$');

            addaxis((CoM_x+lzero)*mag, Energy, '.b', 'LineStyle','none', 'MarkerSize', 16);
            addaxislabel(3, '$Bending\ energy\ (J)$');

            centers(:, 2) = 2048 - centers(:, 2);
            viscircles(centers*mag, radii*mag,'LineStyle','--', 'LineWidth', 0.5, 'Color', 'k'); hold on

            savepath = ['F:\Processing & Results\Actin Filaments in Porous Media\Figures\BendingEnergy_Lee',...
                pathname(56:70), pathname(end-7:end)]; 
            if ~exist(savepath, 'dir')
                mkdir(savepath)
            end

            f=gcf;
            exportgraphics(f,[savepath, filesep, filename(1: end-4), '_E_Lee.png'],'Resolution',100)

            close all
            clearvars CoM_x Energy xy centers radii
        end
    end
end


%%
%%%% Calculate the bending energy and Lee of actin filaments.
%%%% data from Vic_ActinAddInformationSave.m
%%%% saving name format: trajectory_..._batch1_Info.mat

clear; close all; clc;

B = 6.9e-26;  % Bending rigidity

xlsfile = readcell('ForActinPostprocessing.xlsx','Sheet','Sheet1','NumHeaderLines',1);
% This is the file that contains all the information about the later processing (in sheet 1).

NumGroup = size(xlsfile, 1);  % Number of the groups to be calculated.
ExpDate = xlsfile(:, 1);  % The experiment date.
storePath = xlsfile(:, 2);  % Path of the data to be processed.

for no_Group = 23:25

    the_exp_date = yyyymmdd(ExpDate{no_Group, 1});
    thefiles = dir(fullfile(storePath{no_Group},'*.mat'));

    theimgs = dir(['Z:\Experimental Data (EXTRACTED)\Actin Filaments in Porous' ...
        ' Media\',num2str(the_exp_date),'-Actin\AfterAveBGR\*.tif']);

    for file_ind = 1:length(thefiles)

        filename = thefiles(file_ind).name;

        if contains(filename, 'PAsInfoAdded_')

            load(fullfile(thefiles(1).folder, thefiles(file_ind).name));

            pathname = thefiles(1).folder;
            filename = thefiles(file_ind).name

            % find the cooresponding *.tif image
            image_names = struct2cell(theimgs); image_names = image_names(1, :);
            image_ind = find(cellfun(@(x) contains(x, filename(25:end-11)), image_names));

            spl_Ls = xy.arclen_spl(Good_case_frm);
            L_0 = Vic_Get_ave_cutExtreme(spl_Ls, 0.2); % get the filament length
            %     L_0 = Vic_Get_L_histo(spl_Ls); % get the filament length (not good)

            the_fully_stretched_No = find(abs(spl_Ls - L_0) == min(abs(spl_Ls - L_0)));
            base_fiber_no = Good_case_frm(the_fully_stretched_No);
            base_II = imread(fullfile(theimgs(1).folder, theimgs(image_ind).name), xy(1).frame(base_fiber_no)); % load the image

            [intensity_base, ~, ~] = Vic_Fiber_Intensity(xy, the_fully_stretched_No, base_fiber_no, base_II, 30, lzero);

            CoM_x = zeros(size(Good_case_frm,2), 1); Energy = zeros(size(Good_case_frm,2), 1);
            L_ee_norm_belowOne =  zeros(size(Good_case_frm,2), 1); L_ee_norm =  zeros(size(Good_case_frm,2), 1);
            for frm_ind = 1:size(Good_case_frm,2)

                xy_ind = Good_case_frm(frm_ind);% index of the 'good' cases
                II = imread(fullfile(theimgs(1).folder, theimgs(image_ind).name), xy(1).frame(xy_ind)); % load the image

                [~, IntensityAll,R2] = Vic_Fiber_Intensity(xy, frm_ind, xy_ind, II, 30, lzero);

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

                    L_ee_norm(frm_ind) = sqrt((xy.spl{xy_ind}(1,1)-xy.spl{xy_ind}(end,1))^2 ...
                        + (xy.spl{xy_ind}(1,2)-xy.spl{xy_ind}(end,2))^2) / L_0;
                    L_ee_norm_belowOne(frm_ind) = L_ee_norm(frm_ind);
                    if L_ee_norm_belowOne(frm_ind) > 1
                        L_ee_norm_belowOne(frm_ind) = 1;
                    end

                    R2(R2<5) = nan; % remove the too-high curvature part (here, the allowed minimum curvature radius is 0.5um).

                    Energy(frm_ind) = B / 2 * sum((1./R2(2:end)).^2 .* seglen_spl, ...
                        'omitnan') / 1e-7 + fold_no * (B / 2 * (1e6)^2 * pi * 1e-6);
                    % 1e-7 is the scale m/pixel;
                end
            end

            save([pathname, filesep, filename(14: end-4), '_Info.mat'],'centers', ...
                'circleMask', 'CoM_x', 'Energy', 'Good_case_frm', 'L_ee_norm', ...
                'L_ee_norm_belowOne', 'lzero','metric','radii','xy')

            clearvars CoM_x Energy xy centers radii L_ee_norm_belowOne L_ee_norm
        end
    end
end







%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function addaxis(varargin)
%ADDAXIS  adds an axis to the current plot

%  NOTE:  Zhibo modified version (20230207)

%  get current axis
cah = gca;

%  new for R2010a
hzoom = zoom;
set(hzoom,'ActionPostCallback',@addaxis_zoom_post);
set(hzoom,'ActionPreCallback',@addaxis_zoom_pre);

if nargin>=3 && ~ischar(varargin{3})
    yl2 = varargin{3};
    indkeep = setdiff(1:nargin,3);
    [varargintemp{1:nargin-1}] = deal(varargin{indkeep});
    varargin = varargintemp;
end

%  assume existing plot has axes scaled the way you want.
yl = get(cah,'ylim');
cpos = get(cah,'position');
set(cah,'box','off');

%  get userdata of current axis.  this will hold handles to
%  additional axes and the handles to their corresponding plots
%  in the main axis
%  axh = get(cah,'userdata');
axh = getaddaxisdata(cah,'axisdata');

ledge = cpos(1);
if length(axh)>=1
    if length(axh)/2 == round(length(axh)/2)
        rpos = get(axh{end-1}(1),'position');
        redge = rpos(1);
        lpos = get(axh{end}(1),'position');
        ledge = lpos(1);
    else
        rpos = get(axh{end}(1),'position');
        redge = rpos(1);
        if length(axh)>1
        	lpos = get(axh{end-1}(1),'position');
        	ledge = lpos(1);
        end
    end
else
    redge = cpos(3)+cpos(1);
    ledge = cpos(1);
end

totwid = redge-ledge;

%  assume axes are added on right, then left, then right, etc.
numax = length(axh)+1;

%%%%%%%%%%%%% Here change the separation distance %%%%%%%%%
%  parameters setting axis separation
axcompleft=0; 
if numax == 1
    axcompright = 0.0;
else
    axcompright = 0.12;
end

if numax/2 == round(numax/2)
    side = 'left';
    xpos = ledge-axcompleft*totwid;
else
    side = 'right';
    xpos = redge+axcompright*totwid;
end

h_ax = axes('position',[xpos, cpos(2), cpos(3)*.015, cpos(4)]);
%  plot in new axis to get the automatically generated ylimits
hplt = plot(varargin{:});

if ~exist('yl2')
    yl2 = get(h_ax,'ylim');
end


set(h_ax,'yaxislocation',side);
set(h_ax,'color',get(gcf,'color'));
set(h_ax,'box','off');
set(h_ax,'xtick',[]);
set(hplt,'visible','off');

set(h_ax,'ylim',yl2);


%  rescale all y-values
y = varargin{2};

y = (y-yl2(1))./(yl2(2)-yl2(1)).*(yl(2)-yl(1))+yl(1);

varargin{2} = y;
axes(cah)
hplts = aa_splot(varargin{:});
set(gca,'ylim',yl);

%  store the handles in the axis userdata
axh{length(axh)+1} = [h_ax;hplts];
% set(cah,'userdata',axh);
setaddaxisdata(cah,axh,'axisdata');
set(cah,'box','off');

%  set the axis color if a single line was added to the plot
if length(hplts)==1
    set(h_ax,'ycolor',get(hplts,'color'));
end
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



function [intensity_0, IntensityAll, R] = Vic_Fiber_Intensity(xy, frm_ind, xy_ind, II, cut_length, lzero)
% get the intensity along the fiber
% xy -- the struct contains the reconstructed information
% frm_ind -- index of current Good_case_frm (only for lzero)
% xy_ind -- current index of xy
% II -- the corresponding raw image
% cut_length -- the length of normal cut (should be even number)
% lzero -- a useless parameter (but necessary).
%
% IntensityAll -- the integral along the centerline of the filament
% intensity_0 -- the average or common value along the centerline of the filament
% R2 -- Radius of curvature along the fiber

spl = xy(1).spl{xy_ind}; % the x-y coordinates of the spline (in 'plot' coordinate)
spl(:, 1) = spl(:, 1) + lzero(frm_ind); % lzero is very important !!!
spl(:, 2) = size(II, 2) - spl(:, 2) - lzero(frm_ind); % convert the x-y coordinates to 'image' coordinate

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
    %     plot(xy.spl{xy_ind}(:,1)+lzero(frm_ind),size(II, 2)-xy.spl{xy_ind}(:,2)-lzero(frm_ind), 'LineWidth', 2);hold on
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