clear; close all; clc;

magni = 0.1; % the magnification of the objective (um/pixel)
crop_size = 512;
ds = 5; % constant segment length (in px) used for spacing the reference points in the B-spline fitting

% The .mat file where stores your results.
mat_pathname = dir(['Z:\Processing & Results\Actin tracking\' ...
    'Z-scan data\20240223-Z_scan\results\*.mat']);  
rawdata_pathname = 'Z:\Experimental Data (RAW)\Actin tracking\Z-scan data\20240223-Z_scan';
save_pathname_mother = 'Z:\Processing & Results\Actin tracking\Z-scan data\20240223-Z_scan-Reconstruction\';

for ii = 1:length(mat_pathname)

    mat_filename = mat_pathname(ii).name;
    load(fullfile(mat_pathname(1).folder, mat_filename));

    % The corresponding RAW data folder.
    rawdata_foldername = extractBetween(mat_filename,'trajectory_','_batch1');
    % The corresponding RAW tifs
    rawdata_tifs = dir(fullfile(rawdata_pathname,rawdata_foldername{1},'*.tif'));
    img_no_str = cat(1, rawdata_tifs(:).name);
    img_no = str2num(img_no_str(:,6:12)); % extract the image number
    % The corresponding stage position infomation.
    StagePosition_txt = readmatrix(fullfile(rawdata_pathname,rawdata_foldername{1},'StagePosition.txt'),'Range','A:D');
    StagePosition_txt(~ismember(StagePosition_txt(:,1), img_no), :) = []; % remove the invalid rows
    StagePosition_txt = [StagePosition_txt(1, :); StagePosition_txt];
    StagePosition_txt(1) = StagePosition_txt(1) - 1; % add the missed FIRST frame
    if ~isempty(img_no(~ismember(img_no, StagePosition_txt(:,1))))
        missed_frame = find(img_no == img_no(~ismember(img_no, StagePosition_txt(:,1)))); % find the missed frame
    else
        missed_frame = [];
    end
    if length(img_no) < size(StagePosition_txt, 1) + length(missed_frame)
        StagePosition_txt(diff(StagePosition_txt(:,1))==0,:) = []; % remove the duplicate
    end

    Good_case(ismember(Good_case, missed_frame))= []; % no information about the missed frame

    % Only choose the good cases
    for xy_no = 1:length(xy)

        false_ind = ~ismember(xy(xy_no).frame, Good_case);
        % Here, the 'Good_case' means the absolute frame number, not the index of
        % xy.frame (different to the code before which calculation and selection
        % are two processes.

        xy(xy_no).crd(false_ind) = [];
        xy(xy_no).centroid(false_ind) = [];
        xy(xy_no).arclen(false_ind) = [];
        xy(xy_no).seglen(false_ind) = [];
        xy(xy_no).nGoodframe = length(Good_case) - xy(1).nframe + xy(xy_no).nframe;
        xy(xy_no).frame(false_ind) = [];
        xy(xy_no).spl(false_ind) = [];
        xy(xy_no).knots(false_ind) = [];
        xy(xy_no).arclen_spl(false_ind) = [];
        xy(xy_no).seglen_spl(false_ind) = [];

    end

    img_no_Good_case = img_no(Good_case, :);
    StagePosition_Good_case = StagePosition_txt;
    StagePosition_Good_case(~ismember(StagePosition_Good_case(:,1), ...
        img_no_Good_case), :) = [];

    % put information together (good cases only)
    StagePosition_and_frameNo = [StagePosition_Good_case, Good_case'];

    crd_xyz = [];
    for xy_no = 1:length(xy)
        for j = 1:xy(xy_no).nGoodframe

            crd_x = xy(xy_no).crd{j}(:,1); crd_y = xy(xy_no).crd{j}(:,2); 
            if length(crd_x) == 1
                crd_x = [crd_x-2; crd_x; crd_x+2; crd_x];
                crd_y = [crd_y; crd_y-2; crd_y; crd_y+2]; % enlarge the area for single point!            
            end
            crd_y = 1024-crd_y; % 1024 is the height of the image (to image coordinates)

            current_abs_frame = xy(xy_no).frame(j);

            current_crd_z = StagePosition_and_frameNo(StagePosition_and_frameNo(:, 5) == current_abs_frame, 4);
            current_crd_xyz = [crd_x, crd_y, ones(length(crd_x),1)*current_crd_z/magni]; % unit:pixel
            crd_xyz = [crd_xyz; current_crd_xyz]; % Notice the coordinate!
 
        end
    end

    center_xyz = round(mean(crd_xyz, 1));
    crd_xyz_int = round(crd_xyz);

    if min(crd_xyz_int(:,3)) > 0
        the_tmplt = zeros(max(crd_xyz_int(:,2))+10, max(crd_xyz_int(:,1))+10, max(crd_xyz_int(:,3))+10);
        the_tmplt(sub2ind(size(the_tmplt), crd_xyz_int(:,2), crd_xyz_int(:,1), crd_xyz_int(:,3))) = 1;
        % pad the xyz in a 'box'.
        im_blur = imgaussfilt3(the_tmplt,3); % Guassian filter
%         volshow(im_blur);
        im_skel = bwskel(imbinarize(im_blur),'MinBranchLength',10); % skeletonized in 3D
%         volshow(im_skel, 'BackgroundColor','b')
        [y_skel, x_skel, z_skel] = ind2sub(size(the_tmplt), find(im_skel == 1)); % xyz of the skeleton
        skel_xyz = [x_skel, y_skel, z_skel];
    else
        shift_up = abs(min(crd_xyz_int(:,3))) + 10;
        crd_xyz_int(:,3) = crd_xyz_int(:,3) + shift_up; % keep index positive
        the_tmplt = zeros(max(crd_xyz_int(:,2))+10, max(crd_xyz_int(:,1))+10, max(crd_xyz_int(:,3))+10);
        the_tmplt(sub2ind(size(the_tmplt), crd_xyz_int(:,2), crd_xyz_int(:,1), crd_xyz_int(:,3))) = 1;
        % pad the xyz in a 'box'.
        im_blur = imgaussfilt3(the_tmplt,3); % Guassian filter
%         volshow(im_blur);
        im_skel = bwskel(imbinarize(im_blur),'MinBranchLength',10); % skeletonized in 3D
%         volshow(im_skel, 'BackgroundColor','b')
        [y_skel, x_skel, z_skel] = ind2sub(size(the_tmplt), find(im_skel == 1)); % xyz of the skeleton
        z_skel = z_skel - shift_up; % shift back
        skel_xyz = [x_skel, y_skel, z_skel];
    end

    skel_xy_proj = logical(sum(im_skel, 3));
    end_pts = find_skel_ends(skel_xy_proj,'not testing'); % This is in image coordinates
    [~,ipr] = max(sqrt(end_pts(:,1).^2 + end_pts(:,2).^2));
    [~,pos] = ismember(end_pts(ipr,:),[x_skel, y_skel],'rows');



    %%%%%%%%%%%%%%%%%%%%% sort coordinates (start) %%%%%%%%%%%%%%%%%%%%%%%%
    % define the coordinates at the one extremity of the filament
    xsrt = skel_xyz(pos,1);
    ysrt = skel_xyz(pos,2);
    zsrt = skel_xyz(pos,3);

    % rearrange the coordinate matrix by putting the starting point at the beginning of the matrix
    strpnts = [xsrt, ysrt, zsrt];
    skel_xyz(pos, :) = [];
    skel_xyz = [strpnts; skel_xyz];

    % compute the distance matrix between all coordinates
    Dij = squareform(pdist(skel_xyz,'Euclidean')); % D(i,j) is the distance between points i and j
    Dij(Dij==0)=NaN;  % replace 0 by NaN for exluding O from the minima

    numcrd = size(skel_xyz,1); % number of points in the filament
    XYZ = strpnts; % define the matrix where to store the sequential coordinates

    k=1;
    ind=1;  % start appending the coordinates from the starting point
    Dij(:,ind) = NaN; % To avoid re-counting the (first) reference point

    while k < numcrd

        [~,ind] = min(Dij(ind,:)); % point having the smallest distance from the reference point
        Dij(:,ind) = NaN; % replace NaN to avoid counting the same point two times

        XYZ = [XYZ; skel_xyz(ind,:)];  % append the next coordinates sequentially
        k=k+1;

    end
    %%%%%%%%%%%%%%%%%%%%% sort coordinates (end) %%%%%%%%%%%%%%%%%%%%%%%%%%



    %%%%%%%%%%%%%%%%%%%%% spline centerline (start) %%%%%%%%%%%%%%%%%%%%%%%
    [arcl,segl] = arclength(XYZ(:,1),XYZ(:,2),XYZ(:,3),'linear');

    cum = cumsum(segl(:)); % array of the cumulative sum of centerline segments (cum(end)==arclen!!)
    % remainder of the division between cum and ds.
    % Minima are found when the partial sum of the segments is a multiple of ds

    modulo = mod(cum,ds);
    % obtain these local minima and the corresponding array of index. For these indexes,
    % the arclength is an integer multiple of ds. The corresponding  centerline coordinates are thus equispaced

    posmin = find(islocalmin(modulo)==1);
    nend = size(XYZ, 1); % number of coordinates
    if pdist(XYZ(nend-1:nend, :),'Euclidean') < ds
        posmin = [1; posmin+1; nend];
    else
        posmin = [1; posmin+1]; %
    end

    knots = [XYZ(posmin,1), XYZ(posmin,2), XYZ(posmin,3)];
    spline_xyz = BSpline(knots,'order',2,'nint',5)*magni; % unit:um
    spline_xyz_in_pixel = spline_xyz/magni; % unit:pixel; 
    %%%%%%%%%%%%%%%%%%%%% spline centerline (end) %%%%%%%%%%%%%%%%%%%%%%%%%%



    %%%%%%%%%%%%%%%%%%%%% plot to check (start) %%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Color','w','Position',[100 100 800 600]);
    plot3((x_skel-min(x_skel))*magni, (y_skel-min(y_skel))*magni, ...
        (z_skel-min(z_skel))*magni,'r.', 'LineStyle', 'none','MarkerSize',18)
    xlabel('$x\ (\mathrm{\mu m})$','FontSize', 18,'Interpreter', 'latex');
    ylabel('$y\ (\mathrm{\mu m})$','FontSize', 18,'Interpreter', 'latex');
    zlabel('$z\ (\mathrm{\mu m})$','FontSize', 18,'Interpreter', 'latex');
    axis equal
    set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'ZGrid', 'On', 'GridAlpha', 0.5, 'FontSize', 18)

    hold on; 
    plot3(spline_xyz(:,1)-min(spline_xyz(:,1)), spline_xyz(:,2)-min(spline_xyz(:,2)), ...
        spline_xyz(:,3)-min(spline_xyz(:,3)),'k','LineWidth',2)
    %%%%%%%%%%%%%%%%%%%%% plot to check (end) %%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %%%%%%%%%%%%%%%% crop and save -- reconstruction (start) %%%%%%%%%%%%%%%
    for iii = 1:length(StagePosition_and_frameNo)

        image_ind = find(img_no == StagePosition_and_frameNo(iii,1)); % find the cooresponding image
        II = imread(fullfile(rawdata_tifs(1).folder, rawdata_tifs(image_ind).name)); % load the image

        if round(center_xyz(1)-crop_size/2)>0 && round(center_xyz(2)-crop_size/2)>0 && ...
                round(center_xyz(1)+crop_size/2)<size(II,2)+1 && round(center_xyz(2)+crop_size/2)<size(II,1)+1

            Cropped_II = II(round(center_xyz(2)-crop_size/2):round(center_xyz(2)+crop_size/2)-1, ...
                round(center_xyz(1)-crop_size/2):round(center_xyz(1)+crop_size/2)-1);
            if ~isa(Cropped_II,'uint8')
                Cropped_II_uint8 = uint8(double(Cropped_II)/double(max(max(Cropped_II))) * 255);
                % Convert to 8-bits by linearly scaling from min-max to 0-255 (similar to ImageJ)
            else
                Cropped_II_uint8 = Cropped_II;
            end

            if ~exist(save_pathname_mother,"dir")
                mkdir(save_pathname_mother);
            end
            save_tif_name = fullfile(save_pathname_mother, [rawdata_foldername{1}, ...
                '_', rawdata_tifs(image_ind).name]);
            imwrite(Cropped_II_uint8, save_tif_name);
        end

        Z_in_image = StagePosition_and_frameNo(iii,4)/magni; % unit:pixel;

        spline_xyz_shifted = spline_xyz_in_pixel;
        spline_xyz_shifted(:,3) = spline_xyz_shifted(:,3) - Z_in_image; % unit:pixel;
        spline_xyz_shifted(:,1:2) = spline_xyz_shifted(:,1:2) - crop_size/2;

        fileID = fopen([save_tif_name(1:end-4),'_FocusPosition.txt'],'w');
        fprintf(fileID,'%3d %3d %3d\n', spline_xyz_shifted'); % image coordinates!!!
        fclose(fileID);

    end
    %%%%%%%%%%%%%%% crop and save -- reconstruction (end) %%%%%%%%%%%%%%%%%

    
    
    %%%%%%%%%%%%%%%%%% crop and save -- tracking (start) %%%%%%%%%%%%%%%%%
    delta_z = abs(StagePosition_and_frameNo(:,4)/magni - center_xyz(3)); % unit:pixel
    frameNo_Z_center = StagePosition_and_frameNo(delta_z == min(delta_z), 5);
    FocusPosition = [crop_size/2-1, crop_size/2-1, frameNo_Z_center(1)];

    % crop the images and save to a new folder
    for image_ind = 1:length(rawdata_tifs)

        II = imread(fullfile(rawdata_tifs(1).folder, rawdata_tifs(image_ind).name)); % load the image

        if round(center_xyz(1)-crop_size/2)>0 && round(center_xyz(2)-crop_size/2)>0 && ...
            round(center_xyz(1)+crop_size/2)<size(II,2)+1 && round(center_xyz(2)+crop_size/2)<size(II,1)+1

            Cropped_II = II(round(center_xyz(2)-crop_size/2):round(center_xyz(2)+crop_size/2)-1, ...
                round(center_xyz(1)-crop_size/2):round(center_xyz(1)+crop_size/2)-1);
            if ~isa(Cropped_II,'uint8')
                Cropped_II_uint8 = uint8(double(Cropped_II)/double(max(max(Cropped_II))) * 255);
                % Convert to 8-bits by linearly scaling from min-max to 0-255 (similar to ImageJ)
            else
                Cropped_II_uint8 = Cropped_II;
            end

            save_pathname = fullfile(save_pathname_mother, rawdata_foldername{1});
            if ~exist(save_pathname,"dir")
                mkdir(save_pathname);
            end
            imwrite(Cropped_II_uint8, fullfile(save_pathname, rawdata_tifs(image_ind).name));
        end
    end

    copyfile(fullfile(rawdata_pathname,rawdata_foldername{1},'StagePosition.txt'), ...
        fullfile(save_pathname,'StagePosition.txt'));
    copyfile(fullfile(rawdata_pathname,rawdata_foldername{1},'StageCommand.txt'), ...
        fullfile(save_pathname,'StageCommand.txt'));
    copyfile(fullfile(rawdata_pathname,rawdata_foldername{1},'RawImageInfo.txt'), ...
        fullfile(save_pathname,'RawImageInfo.txt'));

    fileID = fopen(fullfile(save_pathname,'FocusPosition.txt'),'w');
    fprintf(fileID,'%3d %3d %3d\n', FocusPosition);
    fclose(fileID);
    %%%%%%%%%%%%%%%%%% crop and save -- tracking (end) %%%%%%%%%%%%%%%%%%%

end



function BS = BSpline(knots,varargin)
%BSPLINE computes the B-spline approximation from a set of coordinates.
% BSPLINE(KNOTS) returns the B-spline interpolation between the reference
% points (knots) whose coordinates are given by the array KNOTS.
% The coordinates of the knots are given vertically, i.e. KNOTS(i,j) gives
% the j-th coordinate of the i-th knot. The knots can be of any dimension.
%
% BSPLINE(KNOTS,'order',n) Uses n -th order approximation (default: n=2)
%
% BSPLINE(KNOTS,'nint',m) gives m points per interval (default: m=10) 
% Zhibo: include the two endpoints.
% 
% If KNOTS is of size [p,q], the result will be of size [(p-1)*(m-1)+1 ,q],
% except if periodicity is requested (see below).
%
% BSPLINE(KNOTS,'periodic',true) use periodic conditions at end knots.
%
ip = inputParser;
addOptional(ip,'order',2)
addOptional(ip,'nint',10)
addOptional(ip,'periodic',false)
parse(ip,varargin{:});

if ip.Results.periodic
    np_rep=ip.Results.order+1;
    knots=[knots(end-np_rep+1:end,:); knots; knots(1:np_rep,:)];
end

p=size(knots,1);
q=size(knots,2);

if p<=2
    BS=knots;
    return
end

n=ip.Results.nint;
n=(n-1)*(p-1)+1;	% Overall number of queried points
y = linspace(0,1,n);
order=min(ip.Results.order+1,p);
Xl = zeros(order,order,q);
t = [zeros(1,order-1),linspace(0,1,p-order+2),ones(1,order-1)]; % node vector £¨= number of knots + number of the orders (not degrees)£©
BS=zeros(n,q);
m_min=1;
m_max=n;
for m = 1:n-1
    t0 = y(m);
    k = find(t0 >= t,1,'last');
    if (k > p)
        BS=BS(1:m-1,:);
        return;
    end
    Xl(:,1,:) = knots(k-order+1:k,:);   
    if k<=order+1
        m_min=max(m,m_min);
    end
    if k>=p-order+2
        m_max=min(m,m_max);
    end
    for i = 2:order
        for j = i:order
            num = t0-t(k-order+j);
            if num == 0
                wt = 0;
            else
                s = t(k+j-i+1)-t(k-order+j);
                wt = num/s;
            end
            Xl(j,i,:) = (1-wt)*Xl(j-1,i-1,:) + wt*Xl(j,i-1,:);
        end
    end
    BS(m,:) = Xl(order,order,:);
end
BS(end,:)=knots(end,:);
if ip.Results.periodic
    BS=BS(m_min:m_max-1,:);
    BS(end,:)=BS(1,:);
end
end
