%%%% Plot the chronophotograoh, bending energy and Lee of actin filaments.
%%%% data from Vic_ActinChronophotograph_BendingE and Vic_ActinChronophotograph_Lee
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

for no_Group = [7 8 13 14 15 16 17 18]

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