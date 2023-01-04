%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Good_case: stored the absolute slice numbers of the calculated *.tif
%             files. Get directly after running the 'Filament-detecting'
%             (update parameters version). 
%             * call in loop j: xy(1).frame(j)
%
%  Good_case_frm: stored the index of the 'Good_case' variable. Get from
%                 the 'draw reconstruction and selection' process.
%                 * call in loop j: xy(1).frame(Good_case_frm(j))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Make a video (CoM fixed) ...
clear; close all; clc;

% The .tif file you just calculated.
[filename, pathname]=uigetfile({['F:\Experimental Data (EXTRACTED)\FSI - ' ...
    'Rigid Fiber &  Individual Obstacle\20221004-SU8_Fibers-Individual_triang' ...
    'ularPillar_uppoint\Inverted\*.tif']}, 'Choose a *.tif to be processed');  % input file
basepath=pathname;
tifname=filename;

% The .mat file where stores your results.
[filename, pathname]=uigetfile({['F:\Processing & Results\FSI - Rigid Fiber ' ...
    '&  Individual Obstacle\20221004-SU8_Fibers-Individual_triangularPillar_' ...
    'uppoint\results\*.mat']}, 'Choose a *.mat to be processed');  % input file
load([pathname, filename]);
close all;
lzero = 0;
xwin = prmt(1).xwin; ywin = prmt(1).ywin;

v = VideoWriter(strcat(filename),'MPEG-4');   % To make a video!
v.FrameRate = 10;  % Frame rate in the video.
v.Quality = 100;

open(v);   % For the video.
for j = 1:size(Good_case_frm,2)

    i = Good_case_frm(j);

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
        round(max(T_C1-xwin/2,1)):round(min(T_C1+xwin/2,size(Foo,2))))*2;

    the_frame = uint8(zeros(ywin+1, xwin+1));  % Put the Zoom_in into this frame so that every figure has the same size to make the video.
    xy_mov = size(the_frame)-size(Zoom_in);  % How much should the plot be shifted.
    if abs(j - 1)  < abs(j - size(Good_case_frm,2))  % To fix the center-of-mass position in the video.
        the_frame(padarray(true(size(Zoom_in)), xy_mov, 'pre')) = Zoom_in;  % insert the Zoom_in matrix into the frame.

        imshow(the_frame, 'InitialMagnification', 200); hold on;
        % Show the B-spline, the idea is to change the original point according
        % to the required window.
        % h = plot(T_spl1,T_spl2,'-','linewidth',6);
        h = plot(T_spl1-max(T_C1-xwin/2,1) + xy_mov(2), T_spl2-max(T_C2-ywin/2,1) + xy_mov(1),'-m','linewidth', 2);
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
        h = plot(T_spl1-max(T_C1-xwin/2,1), T_spl2-max(T_C2-ywin/2,1) + xy_mov(1),'-m','linewidth', 2);
        set(h,'marker','.');
        title(['No.',num2str(j)],'Color','red','FontSize',14);
        hold on
        plot(T_spl1(RR2==min(RR2))-max(T_C1-xwin/2,1),T_spl2(RR2==min(RR2))...
            -max(T_C2-ywin/2,1) + xy_mov(1),'o','markeredgecolor','g','markerfacecolor','g','MarkerSize',5);
        hold on;
        plot(T_C1-max(T_C1-xwin/2,1),T_C2-max(T_C2-ywin/2,1) + xy_mov(1),'*','markeredgecolor','r','MarkerSize',5);
        hold off;
    end

    %     saveas(gcf,'F:\Code\tmpData\tmp','tiffn');
    %     imwrite(imread('F:\Code\tmpData\tmp.tif'), [basepath,'results\Batch_1',tifname], 'writemode', 'append');

    %     plot(lzero+xy(1).spl{i}(:,1),size(Foo,1)*1.5-lzero-xy(1).spl{i}(:,2),'-','linewidth',6)
    %     hold on
    %     plot(lzero+xy(1).crd{i}(:,1),size(Foo,1)*1.5-lzero-xy(1).crd{i}(:,2),'.','markeredgecolor','k','markerfacecolor','w','linewidth',2)
    %     hold on

%     pause(0.001);
    frame = getframe(gcf);
    writeVideo(v,frame);
    close

end
close(v);  % For the video.

% figure('color', 'w'); set(gcf, 'Position', [100 300 1000 500]);
% bgim=imread('D:\Dropbox\tmp\Tri_1 - Copy.tif');  % background image.
% imshow(bgim, []); hold on;

% for j = 1:size(Good_case_frm,2)
% 
%     i = Good_case_frm(j);
%     spl1 = xy(1).spl{i}(:,1); spl2 = xy(1).spl{i}(:,2); % x-y coordinates of the B-spline
%     T_spl1 = lzero+spl1;  T_spl2 = 2048-lzero-spl2;
%     plot(T_spl1,T_spl2)
%     hold on
% 
%     axis equal
%     xlim('auto')
%     ylim('auto')
%     xlabel(' x [ px ] ')
%     ylabel(' y [ px ] ')
% 
% end



%% Make videos for all and wait for selecting good cases ...
clear; close all; clc;

the_folder = dir('F:\Experimental Data (EXTRACTED)\FSI - Rigid Fiber &  Individual Obstacle\20230102-SU8_Fibers-Individual_triangularPillar_uppoint\Inverted\*.tif');
basepath=the_folder(1).folder;
for iiii = 1:length(the_folder)

    tifname=the_folder(iiii).name;
    % The .mat file where stores your results.
    pathname = 'F:\Processing & Results\FSI - Rigid Fiber &  Individual Obstacle\20230102-SU8_Fibers-Individual_triangularPillar_uppoint\results\';
    % if exist([pathname, 'trajectory_', tifname(1:end-4), '_batch1.mat'], 'file')
    if isfile([pathname, 'trajectory_', tifname(1:end-4), '_AABGR_batch1.mat'])
        load([pathname, 'trajectory_', tifname(1:end-4), '_AABGR_batch1.mat']);
        close all;
        lzero = 0;
        xwin = prmt(1).xwin; ywin = prmt(1).ywin;

        v = VideoWriter(strcat(tifname),'MPEG-4');   % To make a video!
        v.FrameRate = 10;  % Frame rate in the video.
        v.Quality = 100;

        open(v);   % For the video.
        for j = 1:size(Good_case,2)

            Foo = imread(fullfile(basepath, tifname),xy(1).frame(j));

            C1 = xy(1).centroid{j}(:,1); C2 = xy(1).centroid{j}(:,2); % Center-of-mass
            T_C1 = lzero+C1;  T_C2 = size(Foo,1)-lzero-C2;
            % Count in the Isero and tranfrom the coordinate because of the different
            % original point between the image and plot
            spl1 = xy(1).spl{j}(:,1); spl2 = xy(1).spl{j}(:,2); % x-y coordinates of the B-spline
            T_spl1 = lzero+spl1;  T_spl2 = size(Foo,1)-lzero-spl2;
            % Same as above

            % calculate the curvature
            X = xy(1).spl{j};
            [L2,R2,K2] = curvature(X);
            RR2 = movmean(R2,ceil(size(X,1)/100));

            % Zoom in (based on the center-of-mass [C1, C2]) and show the image.
            % The size of the ROI depends on the [xwin, ywin]
            fig = figure('Name','filaments','Position', [0 0 100 100]);
            %     imshow(Foo*1000); hold on;
            Zoom_in = Foo(round(max(T_C2-ywin/2,1)):round(min(T_C2+ywin/2,size(Foo,1))),...
                round(max(T_C1-xwin/2,1)):round(min(T_C1+xwin/2,size(Foo,2))))*2;

            the_frame = uint16(zeros(ywin+1, xwin+1));  % Put the Zoom_in into this frame so that every figure has the same size to make the video.
            xy_mov = size(the_frame)-size(Zoom_in);  % How much should the plot be shifted.
            if abs(j - 1)  < abs(j - size(Good_case,2))  % To fix the center-of-mass position in the video.
                the_frame(padarray(true(size(Zoom_in)), xy_mov, 'pre')) = Zoom_in;  % insert the Zoom_in matrix into the frame.

                imshow(the_frame,[], 'InitialMagnification', 300); hold on;
                % Show the B-spline, the idea is to change the original point according
                % to the required window.
                % h = plot(T_spl1,T_spl2,'-','linewidth',6);
                h = plot(T_spl1-max(T_C1-xwin/2,1) + xy_mov(2), T_spl2-max(T_C2-ywin/2,1) + xy_mov(1),'-m','linewidth', 2);
                set(h,'marker','.');
                % title(['No.',num2str(xy.frame(j))],'Color','red','FontSize',14);

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
                imshow(the_frame,[], 'InitialMagnification', 300); hold on;
                h = plot(T_spl1-max(T_C1-xwin/2,1), T_spl2-max(T_C2-ywin/2,1) + xy_mov(1),'-m','linewidth', 2);
                set(h,'marker','.');
                title(['No.',num2str(j)],'Color','red','FontSize',14);
                hold on
                plot(T_spl1(RR2==min(RR2))-max(T_C1-xwin/2,1),T_spl2(RR2==min(RR2))...
                    -max(T_C2-ywin/2,1) + xy_mov(1),'o','markeredgecolor','g','markerfacecolor','g','MarkerSize',5);
                hold on;
                plot(T_C1-max(T_C1-xwin/2,1),T_C2-max(T_C2-ywin/2,1) + xy_mov(1),'*','markeredgecolor','r','MarkerSize',5);
                hold off;
            end

            frame = getframe(gcf);
            writeVideo(v,frame);
            close
        end
        close(v);  % For the video.
    end
end

% % % %%%% select good cases based on the videos ...
% % % clear; clc
% % % 
% % % [filename, pathname]=uigetfile({['D:\Dropbox\PROCESS remotely\Processing & Results' ...
% % %     '\FSI - Rigid Fiber &  Individual Obstacle\20220913-SU8_Fibers-Individual' ...
% % %     '_triangularPillar_uppoint\results\*.mat']}, 'Choose a file to be processed');  % input file
% % % load([pathname, filename]);
% % % Good_case_frm = 1: length(Good_case);
% % % 
% % % % to_be_empty = [13, 15, 46, 47, 48, 49, 51, 52, 57, 58, 59, 60, 61, 62, 63, ... 
% % % %     64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 81, 82, 83, 84, 86, 87, 88, 89, ...
% % % %     90, 91, 92, 93, 94, 95, 96, 100, 103, 104, 105, 106, 107, 108, 109, 110, ...
% % % %     111, 112, 113, 114, 115, 117, 118, 130, 131, 132, 133, 134, 135, 136, ...
% % % %     137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, ...
% % % %     151, 152, 153, 154, 155, 156, 157];
% % % to_be_empty = [];
% % % for i = 1:length(to_be_empty)
% % %     Good_case_frm(to_be_empty(end-i+1)) = [];
% % % end
% % % 
% % % save([pathname, filename],'framelist','Good_case','Good_case_frm','xy', ...
% % %     'ROI','InfoImage','prmt','prcs_img');
