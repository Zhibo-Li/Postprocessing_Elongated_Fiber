%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Good_case: stored the absolute slice numbers of the calculated *.tif
%             files. Get directly after running the 'Filament-detecting'
%             (update parameters version). 
%             * Good_case_frm = find(ismember(xy(1).frame, Good_case));
%             * call in loop j: xy(1).frame(Good_case_frm(j))
%
%  Good_case_frm: stored the index of the 'Good_case' variable. Get from
%                 the 'draw reconstruction and selection' process.
%                 * call in loop j: xy(1).frame(Good_case_frm(j))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Draw reconstruction and selection !!!
% This first part of the code is to select the nice cases after calculating
% the shape based on Francesco's code, and then save these cases into 'Good_case'
% Notice the different 'original points' for the image and coordinate
% system in MATLAB

clear; close all; clc;

% The .tif file you just calculated.
[filename, pathname]=uigetfile({['G:\PhD, PMMH, ESPCI\Experimental Data (EXTRACTED)' ...
    '\20220217-Actin\AfterAveBGR\*.tif']}, 'Choose a file to be processed');  % input file
basepath=pathname;
tifname=filename;

% The .mat file where stores your results.
[filename, pathname]=uigetfile({['G:\PhD, PMMH, ESPCI\Processing\20220217-Actin' ...
    '\results\Group_1\*.mat']}, 'Choose a file to be processed');  % input file
load([pathname, filename]);

lzero = max(lobject,ceil(5*lnoise));
jj = 1;

for i = 1:1:xy.nframe

    % Calculate the length difference between two adjacent filaments, if it's too big (> 20pix) then discard this case.
    % If the length of filament is too small, then discard it, too.
    %     if abs(xy.arclen_spl(i)-xy.arclen_spl(i-1)) >= 20 || xy.arclen_spl(i) <= 30
    %         continue
    %     end

    Foo = imread([basepath, tifname],xy(1).frame(i));

    C1 = xy(1).centroid{i}(:,1); C2 = xy(1).centroid{i}(:,2); % Center-of-mass
    T_C1 = lzero+C1;  T_C2 = size(Foo,1)-lzero-C2;
    % Count in the Isero and tranfrom the coordinate because of the different
    % original point between the image and plot
    spl1 = xy(1).spl{i}(:,1); spl2 = xy(1).spl{i}(:,2); % x-y coordinates of the B-spline
    T_spl1 = lzero+spl1;  T_spl2 = size(Foo,1)-lzero-spl2;
    % Same as above

    % calculate the curvature
    X = xy(1).spl{i};
    [L2,R2,K2] = curvature(X);
    RR2 = movmean(R2,ceil(size(X,1)/100));

    % Zoom in (based on the center-of-mass [C1, C2]) and show the image.
    % The size of the ROI depends on the [xwin, ywin]
    fig = figure('Name','filaments','Position', [0 0 100 100]);
    %     imshow(Foo*1000); hold on;
    Zoom_in = Foo(round(max(T_C2-ywin/2,1)):round(min(T_C2+ywin/2,size(Foo,1))),...
        round(max(T_C1-xwin/2,1)):round(min(T_C1+xwin/2,size(Foo,2))))*1000;
    imshow(Zoom_in, 'InitialMagnification', 200); hold on;

%     xwin = 600;
%     ywin = 300;

    % Show the B-spline, the idea is to change the original point according
    % to the required window.
    h = plot(T_spl1-max(T_C1-xwin/2,1), T_spl2-max(T_C2-ywin/2,1),'-','linewidth',6);
    set(h,'marker','.');
    title(['No.',num2str(xy.frame(i))],'Color','red','FontSize',14); hold on

    plot(T_spl1(RR2==min(RR2))-max(T_C1-xwin/2,1),T_spl2(RR2==min(RR2))...
        -max(T_C2-ywin/2,1),'o','markeredgecolor','g','markerfacecolor','g','MarkerSize',5); hold on;

    plot(T_C1-max(T_C1-xwin/2,1),T_C2-max(T_C2-ywin/2,1),'*','markeredgecolor','r','MarkerSize',5); hold off;

    % This is to draw the relative positions between two neighbouring  figures
    % and check we always follow the same filament.
    if i > 1
        figure('Position', [1500 0 100 100]); imshow(Foo); hold on;
        plot(T_C1_pre, T_C2_pre, 'r*','MarkerSize',3);hold on
        plot(T_C1, T_C2, 'co','MarkerSize',3);hold off
        %         rectangle('Position', [T_C1_pre, T_C2_pre, T_C1-T_C1_pre, T_C2-T_C2_pre], 'EdgeColor', 'b', 'FaceColor', 'r', 'LineWidth', 4); hold off;
    end
    T_C1_pre = T_C1; T_C2_pre = T_C2;


    % Select the good cases manually.
    try
        Inputnum = input('The B-spline fits well?: \n 1 = Yes \n 2 = No \n');
        switch Inputnum
            case 1
                Good_case_frm(jj) = i;
                jj = jj + 1;
                close;
            case 2
                close all;
                continue;
            otherwise
                f = msgbox('Be careful and input again!!');  % To avoid other 'wrong' inputs.
                Inputnum = input('\n \n \n \n Be careful and input again!!! \n The B-spline fits well?: \n 1 = Yes \n 2 = No \n');
                switch Inputnum
                    case 1
                        Good_case_frm(jj) = i;
                        jj = jj + 1;
                        close;
                    case 2
                        close all;
                        continue;
                end
                close all;
        end
        close all;
    catch
        f = msgbox('Be careful and input again!!!'); % To avoid no inputs.
        Inputnum = input('\n \n \n \n Be careful and input again!!! \n The B-spline fits well?: \n 1 = Yes \n 2 = No \n');
        switch Inputnum
            case 1
                Good_case_frm(jj) = i;
                jj = jj + 1;
                close;
            case 2
                close all;
                continue;
        end
        close;
    end
    close all;
end

save([pathname, filename],'thickness','structsensitivity','lnoise','lobject','threshold','ds',...
    'npnts','FilNum','initial_frame',...
    'frame_step','final_frame','framelist','improc','InfoImage',...,
    'sensitivity','MinBranchLength','ROI','missed_frames',...
    'xskip','yskip','xwin','ywin',...
    'N_fil','prcs_img','xy','Good_case_frm')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make video (CoM fixed) (for the 'draw reconstruction and selection' process)!!!
clear; close all; clc;

% The .tif file you just calculated.
[filename, pathname]=uigetfile({['F:\Experimental Data (EXTRACTED)\Actin Filaments' ...
    ' in Porous Media\20210413-Actin\AfterAveBGR\*.tif']}, 'Choose a *.tif to be processed');  % input file
basepath=pathname;
tifname=filename;

% The .mat file where stores your results.
[filename, pathname]=uigetfile({['F:\Processing & Results\Actin Filaments in ' ...
    'Porous Media\20210413-Actin\results\*.mat']}, 'Choose a *.mat to be processed');  % input file
load([pathname, filename]);

lzero = max(lobject,ceil(5*lnoise));

v = VideoWriter(strcat(filename),'MPEG-4');   % To make a video!
v.FrameRate = 24;  % Frame rate in the video.
v.Quality = 100;

open(v);   % For the video.
for j = 1:size(Good_case_frm,2)

    i = Good_case_frm(j);% index of the 'good' cases

    Foo = imread([basepath, tifname],xy(1).frame(i));

    C1 = xy(1).centroid{i}(:,1); C2 = xy(1).centroid{i}(:,2); % Center-of-mass
    T_C1 = lzero+C1;  T_C2 = size(Foo,1)-lzero-C2;
    % Count in the Isero and tranfrom the coordinate because of the different
    % original point between the image and plot
    spl1 = xy(1).spl{i}(:,1); spl2 = xy(1).spl{i}(:,2); % x-y coordinates of the B-spline
    T_spl1 = lzero+spl1;  T_spl2 = size(Foo,1)-lzero-spl2;
    % Same as above

    % calculate the curvature
    X = xy(1).spl{i};
    [L2,R2,K2] = curvature(X);
    RR2 = movmean(R2,ceil(size(X,1)/100));

    % Zoom in (based on the center-of-mass [C1, C2]) and show the image.
    % The size of the ROI depends on the [xwin, ywin]
    fig = figure('Name','filaments','Position', [0 0 100 100]);
    %     imshow(Foo*1000); hold on;
    Zoom_in = Foo(round(max(T_C2-ywin/2,1)):round(min(T_C2+ywin/2,size(Foo,1))),...
        round(max(T_C1-xwin/2,1)):round(min(T_C1+xwin/2,size(Foo,2))))*1000;

    the_frame = uint16(zeros(ywin+1, xwin+1));  % Put the Zoom_in into this frame so that every figure has the same size to make the video.
    xy_mov = size(the_frame)-size(Zoom_in);  % How much should the plot be shifted.
    if abs(j - 1)  < abs(j - size(Good_case_frm,2))  % To fix the center-of-mass position in the video.
        the_frame(padarray(true(size(Zoom_in)), xy_mov, 'pre')) = Zoom_in;  % insert the Zoom_in matrix into the frame.

        imshow(the_frame, 'InitialMagnification', 200); hold on;
        % Show the B-spline, the idea is to change the original point according
        % to the required window.
        % h = plot(T_spl1,T_spl2,'-','linewidth',6);
        h = plot(T_spl1-max(T_C1-xwin/2,1) + xy_mov(2), T_spl2-max(T_C2-ywin/2,1) + xy_mov(1),'-','linewidth',6);
        set(h,'marker','.');
        % title(['No.',num2str(xy.frame(i))],'Color','red','FontSize',14);

        title(['No.',num2str(j)],'Color','red','FontSize',14);
        hold on
        %     quiver(lzero+X(:,1),size(Foo,1)-lzero-X(:,2),K2(:,1),-K2(:,2));  % WHY it is '-K2(:,2)'??
        %     hold on;

        plot(T_spl1(RR2==min(RR2))-max(T_C1-xwin/2,1) + xy_mov(2),T_spl2(RR2==min(RR2))...
            -max(T_C2-ywin/2,1) + xy_mov(1),'o','markeredgecolor','g','markerfacecolor','g','MarkerSize',5);
        hold on;

        plot(T_C1-max(T_C1-xwin/2,1) + xy_mov(2),T_C2-max(T_C2-ywin/2,1) + xy_mov(1),'*','markeredgecolor','r','MarkerSize',5);
        hold off;

    else
        the_frame(padarray(true(size(Zoom_in)), xy_mov, 'post')) = Zoom_in;  % insert the Zoom_in matrix into the frame.
        imshow(the_frame, 'InitialMagnification', 200); hold on;
        h = plot(T_spl1-max(T_C1-xwin/2,1), T_spl2-max(T_C2-ywin/2,1) + xy_mov(1),'-','linewidth',6);
        set(h,'marker','.');
        title(['No.',num2str(j)],'Color','red','FontSize',14);
        hold on
        plot(T_spl1(RR2==min(RR2))-max(T_C1-xwin/2,1),T_spl2(RR2==min(RR2))...
            -max(T_C2-ywin/2,1) + xy_mov(1),'o','markeredgecolor','g','markerfacecolor','g','MarkerSize',5);
        hold on;
        plot(T_C1-max(T_C1-xwin/2,1),T_C2-max(T_C2-ywin/2,1) + xy_mov(1),'*','markeredgecolor','r','MarkerSize',5);
        hold off;
    end

    pause(0.001);
    frame = getframe(gcf);
    writeVideo(v,frame);
    close

end
close(v);  % For the video.

for j = 1:size(Good_case_frm,2)

    i = Good_case_frm(j); % index of the 'good' cases
    spl1 = xy(1).spl{i}(:,1); spl2 = xy(1).spl{i}(:,2); % x-y coordinates of the B-spline
    T_spl1 = lzero+spl1;  T_spl2 = 2048-lzero-spl2;
    plot(T_spl1,T_spl2)
    hold on

    axis equal
    xlim('auto')
    ylim('auto')
    xlabel(' x [ px ] ')
    ylabel(' y [ px ] ')

end

% f=gcf;
% exportgraphics(f,[pathname, filename, 'trajectory.png'],'Resolution',200)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make video (CoM fixed) (the update parameters version of 'Filament-detecting' )!
clear; close all; clc;

% The .tif file you just calculated.
[filename, pathname]=uigetfile({['F:\Experimental Data (EXTRACTED)\Actin Filaments' ...
    ' in Porous Media\20210914-Actin\AfterAveBGR\*.tif']}, 'Choose a *.tif to be processed');  % input file
basepath=pathname;
tifname=filename;

% The .mat file where stores your results.
[filename, pathname]=uigetfile({['F:\Processing & Results\Actin Filaments in ' ...
    'Porous Media\20210914-Actin\results\*.mat']}, 'Choose a *.mat to be processed');  % input file
load([pathname, filename]);

v = VideoWriter(strcat(filename),'MPEG-4');   % To make a video!
v.FrameRate = 24;  % Frame rate in the video.
v.Quality = 100;
open(v);   % For the video.

xwin_max = 0; ywin_max = 0;
for ii = 1: length(prmt)-1
    xwin_max = max([prmt(ii+1).xwin, prmt(ii).xwin, xwin_max]);
    ywin_max = max([prmt(ii+1).ywin, prmt(ii).ywin, ywin_max]);
end
ywin_max = round(ywin_max/2)*2; xwin_max = round(xwin_max/2)*2; 
% The frame size in the video and make sure they are even numbers because 
% will need half of the size. 

false_ind = ~ismember(xy(1).frame, Good_case);
% Here, the 'Good_case' means the absolute frame number, not the index of
% xy.frame (different to the code before which calculation and selection
% are two processes.

xy.crd(false_ind) = [];
xy.centroid(false_ind) = [];
xy.arclen(false_ind) = [];
xy.seglen(false_ind) = [];
xy.nframe = length(Good_case);
xy.frame(false_ind) = [];
xy.spl(false_ind) = [];
xy.knots(false_ind) = [];
xy.arclen_spl(false_ind) = [];
xy.seglen_spl(false_ind) = [];

for j = 1:xy.nframe

    % The real parameters used to do the calculation because the first
    % processing frame might not be 1. Here is mainly for xwin and ywin
    % below.
    prmt_i = xy(1).frame(j) - xy(1).frame(1) + 1;
    Foo = imread([basepath, tifname],xy(1).frame(j));

    % tranfrom the coordinate because of the different original point
    % between the image and plot.
    C1 = xy(1).centroid{j}(:,1); C2 = xy(1).centroid{j}(:,2); % Center-of-mass
    T_C1 = C1;  T_C2 = size(Foo,1)-C2;
    spl1 = xy(1).spl{j}(:,1); spl2 = xy(1).spl{j}(:,2); % x-y coordinates of the B-spline
    T_spl1 = spl1;  T_spl2 = size(Foo,1)-spl2;

    % calculate the curvature
    X = xy(1).spl{j};
    [L2,R2,K2] = curvature(X);
    RR2 = movmean(R2,ceil(size(X,1)/100));

    % Zoom in (based on the center-of-mass [C1, C2]) and show the image.
    % The size of the ROI depends on the [xwin, ywin]
    fig = figure('Name','filaments','Position', [0 0 100 100]);
    Zoom_in = Foo(max(T_C2-ywin_max/2,1):min(T_C2+ywin_max/2,size(Foo,1)),...
        max(T_C1-xwin_max/2,1):min(T_C1+xwin_max/2,size(Foo,2)))*1000;

    % How much should the plot be shifted (ensure the center-of-mass is
    % fixed). It's equal to "xy_mov = size(the_frame)-size(Zoom_in)"
    xy_mov = [ywin_max+1, xwin_max+1]-size(Zoom_in);
    
    if abs(j - 1)  < abs(j - size(Good_case,2))  % To fix the center-of-mass position in the video.
        the_frame = uint16(zeros(ywin_max+1, xwin_max+1));
        % Will put the Zoom_in into this frame so that every figure has the same
        % size to make the video.
        the_frame(padarray(true(size(Zoom_in)), xy_mov, 'pre')) = Zoom_in;  % insert the Zoom_in matrix into the frame.

        imshow(the_frame, 'InitialMagnification', 200); hold on;
        % Show the B-spline, the idea is to change the original point according
        % to the required window.
        % h = plot(T_spl1,T_spl2,'-','linewidth',6);
        h = plot(T_spl1-max(T_C1-xwin_max/2,1) + xy_mov(2), T_spl2-max(T_C2-ywin_max/2,1) + xy_mov(1),'-','linewidth',0.1);
        set(h,'marker','.');
        % title(['No.',num2str(xy.frame(j))],'Color','red','FontSize',14);

        title(['No.',num2str(j)],'Color','red','FontSize',14);
        hold on
        %     quiver(lzero+X(:,1),size(Foo,1)-lzero-X(:,2),K2(:,1),-K2(:,2));  % WHY it is '-K2(:,2)'??
        %     hold on;

        plot(T_spl1(RR2==min(RR2))-max(T_C1-xwin_max/2,1) + xy_mov(2),T_spl2(RR2==min(RR2))...
            -max(T_C2-ywin_max/2,1) + xy_mov(1),'o','markeredgecolor','g','markerfacecolor','g','MarkerSize',5);
        hold on;

        plot(T_C1-max(T_C1-xwin_max/2,1) + xy_mov(2),T_C2-max(T_C2-ywin_max/2,1) + xy_mov(1),'*','markeredgecolor','r','MarkerSize',5);
        hold off;

    else
        the_frame = uint16(zeros(ywin_max+1, xwin_max+1));
        % Will put the Zoom_in into this frame so that every figure has the same
        % size to make the video.
        the_frame(padarray(true(size(Zoom_in)), xy_mov, 'post')) = Zoom_in;  % insert the Zoom_in matrix into the frame.
        imshow(the_frame, 'InitialMagnification', 200); hold on;
        h = plot(T_spl1-max(T_C1-xwin_max/2,1), T_spl2-max(T_C2-ywin_max/2,1) + xy_mov(1),'-','linewidth',0.1);
        set(h,'marker','.');
        title(['No.',num2str(j)],'Color','red','FontSize',14);
        hold on
        plot(T_spl1(RR2==min(RR2))-max(T_C1-xwin_max/2,1),T_spl2(RR2==min(RR2))...
            -max(T_C2-ywin_max/2,1) + xy_mov(1),'o','markeredgecolor','g','markerfacecolor','g','MarkerSize',5);
        hold on;
        plot(T_C1-max(T_C1-xwin_max/2,1),T_C2-max(T_C2-ywin_max/2,1) + xy_mov(1),'*','markeredgecolor','r','MarkerSize',5);
        hold off;
    end

    pause(0.001);
    frame = getframe(gcf);
    writeVideo(v,frame);
    close

end
close(v);  % For the video.

figure('color', 'w'); set(gcf, 'Position', [100 300 1000 500]);
Good_case_frm = find(ismember(xy(1).frame, Good_case));
for j = 1:size(Good_case_frm,2)

    i = Good_case_frm(j);% index of the 'good' cases

    plot(xy.spl{i}(:,1),xy.spl{i}(:,2))
    hold on

    axis equal
    xlim('auto')
    ylim('auto')
    xlabel(' x [ px ] ')
    ylabel(' y [ px ] ')

end
