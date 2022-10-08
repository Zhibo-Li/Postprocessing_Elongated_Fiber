% %% Draw and select the good cases!!!
% This first part of the code is to select the nice cases after calculating
% the shape based on Francesco's code, and then save these cases into 'Good_case'
% Notice the different 'original points' for the image and coordinate
% system in MATLAB

clear; close all; clc;

% The .tif file you just calculated.
[filename, pathname]=uigetfile({['Z:\Experimental Data (EXTRACTED)\FSI - ' ...
    'Rigid Fiber &  Individual Obstacle\20221005-SU8_Fibers-Individual_triang' ...
    'ularPillar_uppoint\Inverted\*.tif']}, 'Choose a file to be processed');  % input file
basepath=pathname;
tifname=filename;

% The .mat file where stores your results.
[filename, pathname]=uigetfile({['Z:\Processing & Results\FSI - Rigid Fiber ' ...
    '&  Individual Obstacle\20221005-SU8_Fibers-Individual_triangularPillar_' ...
    'uppoint\results\*.mat']}, ...
    'Choose a file to be processed');  % input file
load([pathname, filename]);
% clearvars Good_case_frm

% lnoise = prmt(1).lnoise; lobject = prmt(1).lobject;
xwin = prmt(1).xwin; ywin = prmt(1).ywin;
close all;
lzero = 0; %max(lobject,ceil(5*lnoise));
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
        round(max(T_C1-xwin/2,1)):round(min(T_C1+xwin/2,size(Foo,2))))*2;
    imshow(Zoom_in, 'InitialMagnification', 200); hold on;

    % Show the B-spline, the idea is to change the original point according
    % to the required window.
    %     h = plot(T_spl1,T_spl2,'-','linewidth',6);
    h = plot(T_spl1-max(T_C1-xwin/2,1), T_spl2-max(T_C2-ywin/2,1),'-m','linewidth', 2);
    set(h,'marker','.');
    title(['No.',num2str(xy.frame(i))],'Color','red','FontSize',14);
    hold on

    %     quiver(lzero+X(:,1),size(Foo,1)-lzero-X(:,2),K2(:,1),-K2(:,2));  % WHY it is '-K2(:,2)'??
    %     hold on;

    plot(T_spl1(RR2==min(RR2))-max(T_C1-xwin/2,1),T_spl2(RR2==min(RR2))...
        -max(T_C2-ywin/2,1),'o','markeredgecolor','g','markerfacecolor','g','MarkerSize',5);
    hold on;

    plot(T_C1-max(T_C1-xwin/2,1),T_C2-max(T_C2-ywin/2,1),'*','markeredgecolor','r','MarkerSize',5);
    hold off;

    % This is to draw the relative positions between two neighbouring  figures
    % and check we always follow the same filament.
    if i > 1
        figure('Position', [1200 0 100 100]); imshow(Foo); hold on;
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


    %     saveas(gcf,'F:\Code\tmpData\tmp','tiffn');
    %     imwrite(imread('F:\Code\tmpData\tmp.tif'), [basepath,'results\Batch_1',tifname], 'writemode', 'append');

    %     plot(lzero+xy(1).spl{i}(:,1),size(Foo,1)*1.5-lzero-xy(1).spl{i}(:,2),'-','linewidth',6)
    %     hold on
    %     plot(lzero+xy(1).crd{i}(:,1),size(Foo,1)*1.5-lzero-xy(1).crd{i}(:,2),'.','markeredgecolor','k','markerfacecolor','w','linewidth',2)
    %     hold on

    close all;

end

notes = 'The interval of "i (Good_case_frm)" is 1 when selecting the good cases (2022-09-20).';

save([pathname, filename],'framelist','Good_case','Good_case_frm','xy', ...
    'ROI','InfoImage','prmt','prcs_img', 'notes');



%% Draw good cases!!!
% % % clear; close all; clc;
% % % 
% % % % The .tif file you just calculated.
% % % [filename, pathname]=uigetfile({['Z:\Experimental Data (EXTRACTED)\FSI - ' ...
% % %     'Rigid Fiber &  Individual Obstacle\20221004-SU8_Fibers-Individual_triang' ...
% % %     'ularPillar_uppoint\Inverted\*.tif']}, 'Choose a file to be processed');  % input file
% % % basepath=pathname;
% % % tifname=filename;
% % % 
% % % % The .mat file where stores your results.
% % % [filename, pathname]=uigetfile({['Z:\Processing & Results\FSI - Rigid Fiber ' ...
% % %     '&  Individual Obstacle\20221004-SU8_Fibers-Individual_triangularPillar_' ...
% % %     'uppoint\results\*.mat']}, 'Choose a file to be processed');  % input file
% % % load([pathname, filename]);
% % % close all;
% % % lzero = 0;
% % % xwin = prmt(1).xwin; ywin = prmt(1).ywin;
% % % 
% % % v = VideoWriter(strcat(filename),'MPEG-4');   % To make a video!
% % % v.FrameRate = 10;  % Frame rate in the video.
% % % v.Quality = 100;
% % % 
% % % open(v);   % For the video.
% % % for j = 1:size(Good_case_frm,2)
% % % 
% % %     i = Good_case_frm(j);% index of the 'good' cases
% % % 
% % %     Foo = imread([basepath, tifname],xy(1).frame(i));
% % % 
% % %     C1 = xy(1).centroid{i}(:,1); C2 = xy(1).centroid{i}(:,2); % Center-of-mass
% % %     T_C1 = lzero+C1;  T_C2 = size(Foo,1)-lzero-C2;
% % %     % Count in the Isero and tranfrom the coordinate because of the different
% % %     % original point between the image and plot
% % %     spl1 = xy(1).spl{i}(:,1); spl2 = xy(1).spl{i}(:,2); % x-y coordinates of the B-spline
% % %     T_spl1 = lzero+spl1;  T_spl2 = size(Foo,1)-lzero-spl2;
% % %     % Same as above
% % % 
% % %     % calculate the curvature
% % %     X = xy(1).spl{i};
% % %     [L2,R2,K2] = curvature(X);
% % %     RR2 = movmean(R2,ceil(size(X,1)/100));
% % % 
% % %     % Zoom in (based on the center-of-mass [C1, C2]) and show the image.
% % %     % The size of the ROI depends on the [xwin, ywin]
% % %     fig = figure('Name','filaments','Position', [0 0 100 100]);
% % %     %     imshow(Foo*1000); hold on;
% % %     Zoom_in = Foo(round(max(T_C2-ywin/2,1)):round(min(T_C2+ywin/2,size(Foo,1))),...
% % %         round(max(T_C1-xwin/2,1)):round(min(T_C1+xwin/2,size(Foo,2))))*2;
% % % 
% % %     the_frame = uint8(zeros(ywin+1, xwin+1));  % Put the Zoom_in into this frame so that every figure has the same size to make the video.
% % %     xy_mov = size(the_frame)-size(Zoom_in);  % How much should the plot be shifted.
% % %     if abs(j - 1)  < abs(j - size(Good_case_frm,2))  % To fix the center-of-mass position in the video.
% % %         the_frame(padarray(true(size(Zoom_in)), xy_mov, 'pre')) = Zoom_in;  % insert the Zoom_in matrix into the frame.
% % % 
% % %         imshow(the_frame, 'InitialMagnification', 200); hold on;
% % %         % Show the B-spline, the idea is to change the original point according
% % %         % to the required window.
% % %         % h = plot(T_spl1,T_spl2,'-','linewidth',6);
% % %         h = plot(T_spl1-max(T_C1-xwin/2,1) + xy_mov(2), T_spl2-max(T_C2-ywin/2,1) + xy_mov(1),'-m','linewidth', 2);
% % %         set(h,'marker','.');
% % %         % title(['No.',num2str(xy.frame(i))],'Color','red','FontSize',14);
% % % 
% % %         title(['No.',num2str(j)],'Color','red','FontSize',14);
% % %         hold on
% % %         %     quiver(lzero+X(:,1),size(Foo,1)-lzero-X(:,2),K2(:,1),-K2(:,2));  % WHY it is '-K2(:,2)'??
% % %         %     hold on;
% % % 
% % %         plot(T_spl1(RR2==min(RR2))-max(T_C1-xwin/2,1) + xy_mov(2),T_spl2(RR2==min(RR2))...
% % %             -max(T_C2-ywin/2,1) + xy_mov(1),'o','markeredgecolor','g','markerfacecolor','g','MarkerSize',5);
% % %         hold on;
% % % 
% % %         plot(T_C1-max(T_C1-xwin/2,1) + xy_mov(2),T_C2-max(T_C2-ywin/2,1) + xy_mov(1),'*','markeredgecolor','r','MarkerSize',5);
% % %         hold off;
% % % 
% % %     else
% % %         the_frame(padarray(true(size(Zoom_in)), xy_mov, 'post')) = Zoom_in;  % insert the Zoom_in matrix into the frame.
% % %         imshow(the_frame, 'InitialMagnification', 200); hold on;
% % %         h = plot(T_spl1-max(T_C1-xwin/2,1), T_spl2-max(T_C2-ywin/2,1) + xy_mov(1),'-m','linewidth', 2);
% % %         set(h,'marker','.');
% % %         title(['No.',num2str(j)],'Color','red','FontSize',14);
% % %         hold on
% % %         plot(T_spl1(RR2==min(RR2))-max(T_C1-xwin/2,1),T_spl2(RR2==min(RR2))...
% % %             -max(T_C2-ywin/2,1) + xy_mov(1),'o','markeredgecolor','g','markerfacecolor','g','MarkerSize',5);
% % %         hold on;
% % %         plot(T_C1-max(T_C1-xwin/2,1),T_C2-max(T_C2-ywin/2,1) + xy_mov(1),'*','markeredgecolor','r','MarkerSize',5);
% % %         hold off;
% % %     end
% % % 
% % %     %     saveas(gcf,'F:\Code\tmpData\tmp','tiffn');
% % %     %     imwrite(imread('F:\Code\tmpData\tmp.tif'), [basepath,'results\Batch_1',tifname], 'writemode', 'append');
% % % 
% % %     %     plot(lzero+xy(1).spl{i}(:,1),size(Foo,1)*1.5-lzero-xy(1).spl{i}(:,2),'-','linewidth',6)
% % %     %     hold on
% % %     %     plot(lzero+xy(1).crd{i}(:,1),size(Foo,1)*1.5-lzero-xy(1).crd{i}(:,2),'.','markeredgecolor','k','markerfacecolor','w','linewidth',2)
% % %     %     hold on
% % % 
% % % %     pause(0.001);
% % %     frame = getframe(gcf);
% % %     writeVideo(v,frame);
% % %     close
% % % 
% % % end
% % % close(v);  % For the video.
% % % 
% % % % figure('color', 'w'); set(gcf, 'Position', [100 300 1000 500]);
% % % % bgim=imread('D:\Dropbox\tmp\Tri_1 - Copy.tif');  % background image.
% % % % imshow(bgim, []); hold on;
% % % 
% % % % for j = 1:size(Good_case_frm,2)
% % % % 
% % % %     i = Good_case_frm(j); % index of the 'good' cases
% % % %     spl1 = xy(1).spl{i}(:,1); spl2 = xy(1).spl{i}(:,2); % x-y coordinates of the B-spline
% % % %     T_spl1 = lzero+spl1;  T_spl2 = 2048-lzero-spl2;
% % % %     plot(T_spl1,T_spl2)
% % % %     hold on
% % % % 
% % % %     axis equal
% % % %     xlim('auto')
% % % %     ylim('auto')
% % % %     xlabel(' x [ px ] ')
% % % %     ylabel(' y [ px ] ')
% % % % 
% % % % end
