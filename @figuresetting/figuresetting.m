classdef figuresetting

    properties (SetAccess = public)

        units                %  preferred unit for the width and height of the figure
        %                      e.g. 'inches', 'centimeters', 'pixels', 'points', 'characters', 'normalized'
        %                      Default is 'centimeters'
        %
        width                % width of the figure in units defined by 'units'
        %                      Default is 14 centimeters

        height              % height of the figure in units defined by 'units'
        %                      Specifying only one dimension sets the other dimension
        %                      to preserve the figure's default aspect ratio.

        font                  %The font name for all the texts on the figure, including labels, title, legend, colorbar, etc.
        %                      Default is 'Times New Roman'

        fontsize             %The font size for all the texts on the figure, including labels, title, legend, colorbar, etc.
        %                      Default is 14pt

        border               %Thin white border around the figure (compatible with export_fig -nocrop)
        %                      'on', 'off'                      Default is 'off'

        minortickX       % Display minor ticks on the X axis
        %                      'on', 'off'                      Default is 'off'

        minortickY        % Display minor ticks on the Y axis
        %                      'on', 'off'                      Default is 'off'

        linewidth        % Sett the width of the border line
        %                      'on', 'off'                      Default is 'off'
    end

    %% Constructor

    methods

        function obj = figuresetting(units,width,height,font,fontsize,border,linewidth,minortickX,minortickY)
            obj.units = units;
            obj.width = width;
            obj.height = height;
            obj.font = font;
            obj.fontsize=fontsize;
            obj.border=border;
            obj.minortickX=minortickX;
            obj.minortickY=minortickY;
            obj.linewidth=linewidth;
        end
    end


    %% Private methods

    methods

        function obj = figure(obj,titre)

            if nargin>1
                s=figure('name',titre);
            elseif nargin==1
                s=figure;
            end

            bgcolor='w';
            % set the background color
            set(s, 'color',bgcolor);
            % set the font and font size
            set(s, 'DefaultTextFontSize', obj.fontsize);
            set(s, 'DefaultAxesFontSize', obj.fontsize);
            set(s, 'DefaultAxesFontName', obj.font);
            set(s, 'DefaultTextFontName', obj.font);

            %%%%%%%%%%% set the figure size

            % set the root unit
            set(0,'Units',obj.units);
            % get the screen size
            scrsz = get(0,'ScreenSize');
            % set the figure unit
            set(s,'Units',obj.units);
            % get the figure's position
            pos = get(s, 'Position');
            old_pos=pos;

            %set the width and height of the figure
            pos(3)=obj.width;
            pos(4)=obj.height;

            %make sure the figure stays in the middle of the screen
            diff=old_pos-pos;

            if diff(3)<0
                pos(1)=old_pos(1)+diff(3)/2;
                if pos(1)<0
                    pos(1)=0;
                end
            end
            if diff(4)<0
                pos(2)=old_pos(2)+diff(4);
                if pos(2)<0
                    pos(2)=0;
                end
            end

            %warning if the given width (or height) is greater than the screen size
            if pos(3)>scrsz(3)
                warning(['Maximum width (screen width) is reached! width=' num2str(scrsz(3)) ' ' obj.units]);
            end

            if pos(4)>scrsz(4)
                warning(['Maximum height (screen height) is reached! height=' num2str(scrsz(4)) ' ' obj.units]);
            end

            %apply the width, height, and position to the figure
            set(s, 'Position', pos);
            if strcmpi(obj.border, 'off')
                set(s,'DefaultAxesLooseInset',[0,0,0,0]);
            end


        end
    end

    %%
    methods

        function obj = axes(obj,ScaleAxisX,limitx,ScaleAxisY,limity,LabelX,LabelY,FontLabel)%,inposition)
            ax = gca;
            set(ax, 'XMinorTick',obj.minortickX,'YMinorTick',obj.minortickY)
            set(ax,'TickLength',[0.020 0.010],'tickdir','in')
%             set(ax,'Position',inposition)
            set(ax,'XScale',ScaleAxisX,'YScale',ScaleAxisY)
            set(ax,'linewidth',obj.linewidth);
            xlabel(LabelX,'fontsize',FontLabel)
            ylabel(LabelY,'fontsize',FontLabel)
            xlim(limitx)
            ylim(limity)
            box on
        end

        function obj = axes_ticks(obj,set_xticks,set_yticks)
            xticks(set_xticks); 
            yticks(set_yticks);
        end

        function obj = interp_font(obj,interp)
            ax = gca;
            lbx = get(ax,'xlabel');
            lby = get(ax,'ylabel');
            set(lbx,'Interpreter',interp);
            set(lby,'Interpreter',interp);
            % % legobj = findobj('type','legend');
            % % legobj.Interpreter = interp;
        end

        function obj = axesyy(obj,currAx,ScaleAxisX,limitx,ScaleAxisY1,limity1,...
                ScaleAxisY2,limity2,LabelX,LabelY1,LabelY2,FontLabel,inposition)

            set(gca,'Position',inposition)

            % axis 1=right
            set(currAx(1),'yColor','k')
            xlabel(currAx(1),LabelX,'fontsize',FontLabel)
            xlim(currAx(1),limitx)
            set(currAx(1), 'XMinorTick',obj.minortickX,'YMinorTick',obj.minortickY)
            set(currAx(1),'TickLength',[0.020 0.010],'tickdir','out')
            set(currAx(1),'xscale',ScaleAxisX,'yscale',ScaleAxisY1)
            set(currAx(1),'linewidth',obj.linewidth);
            set(currAx(1),'YAxisLocation','left','XAxisLocation','bottom');
            ylabel(currAx(1),LabelY1,'fontsize',FontLabel)
            ylim(currAx(1),limity1)
            box off
            % axis 2=left
            set(currAx(2),'yColor','k')
            %xlim(currAx(2),limitx)
            set(currAx(2), 'XMinorTick',obj.minortickX,'YMinorTick',obj.minortickY)
            set(currAx(2),'TickLength',[0.020 0.010],'tickdir','out')
            set(currAx(2),'xscale',ScaleAxisX,'yscale',ScaleAxisY2)
            set(currAx(2),'linewidth',obj.linewidth);
            set(currAx(2),'YAxisLocation','right','XAxisLocation','top','Xticklabel',[]);
            ylabel(currAx(2),LabelY2,'fontsize',FontLabel)
            ylim(currAx(2),limity2)

        end


    end


    %%
    methods

        function obj = legend(obj,x,location,legendfont,orientation)
            if nargin==1
                l=legend([]);
                set(l,'location','northeast','FontSize',obj.fontsize,'Interpreter','latex')
            elseif nargin==2
                l=legend(x);
                set(l,'location','northeast','FontSize',obj.fontsize,'Interpreter','latex')
            elseif nargin==3
                l=legend(x);
                set(l,'location',location,'FontSize',obj.fontsize,'Interpreter','latex')
            elseif nargin==4
                l=legend(x);
                set(l,'location',location,'FontSize',legendfont,'Interpreter','latex')
            elseif nargin==5
                l=legend(x);
                set(l,'location',location,'FontSize',legendfont,'orientation',orientation,'Interpreter','latex')
            end
        end

    end

end

