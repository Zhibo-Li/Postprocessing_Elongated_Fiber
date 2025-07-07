close all;
clc;
commandwindow;
addpath('D:\Dropbox\GitHub\Postprocessing_Elongated_Fiber') % Change this path to the location of the set_plot.m file (just to have nicer LaTeX plots)

boolLoadVelocity = 0;        % If you need to load the velocity field from converged_flow.vtk (It takes some time, so once it's loaded, you can turn this boolean to 0 before running the script again)

if boolLoadVelocity
    clear;
    boolLoadVelocity = 1;
end

boolPlotVelocities = 0;      % If you need to plot the velocity beads and flow vectors
boolPlotGammaProfiles = 0;   % If you need to plot gamma profiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ THE VELOCITY FIELD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if boolLoadVelocity

    % Read the velocity field
    
    pattern_X = '^X_COORDINATES\s+\d*\s+float$';
    pattern_Y = '^Y_COORDINATES\s+\d*\s+float$';
    pattern_Z = '^Z_COORDINATES\s+\d*\s+float$';
    pattern_velocity = 'VECTORS Velocity float';
        
    filename = 'E:\Experimental Data (RAW)\Simulation\data\converged_flow.vtk'; % Change this accordingly (Zhibo).
    
    fid = fopen(filename,'r');
    
    line = fgetl(fid);
    
    % Read X
    disp('Reading X...')
    while numel(regexp(line,pattern_X)) < 1
        line = fgetl(fid);
    end
    Nx = str2double(regexp(line,'\d*','match'));
    X_grid = zeros(Nx,1);
    for i = 1:Nx
        line = fgetl(fid);
        X_grid(i) = str2double(line);
    end
    
    % Read Y
    disp('Reading Y...')
    while numel(regexp(line,pattern_Y)) < 1
        line = fgetl(fid);
    end
    Ny = str2double(regexp(line,'\d*','match'));
    Y_grid = zeros(Ny,1);
    for i = 1:Ny
        line = fgetl(fid);
        Y_grid(i) = str2double(line);
    end
    
    % Read Z
    disp('Reading Z...')
    while numel(regexp(line,pattern_Z)) < 1
        line = fgetl(fid);
    end
    Nz = str2double(regexp(line,'\d*','match'));
    Z_grid = zeros(Nz,1);
    for i = 1:Nz
        line = fgetl(fid);
        Z_grid(i) = str2double(line);
    end
    
    Np = Nx*Ny*Nz;
    frewind(fid);
    
    % Read the velocity
    disp('Reading velocity...')
    
    % Find first line of the velocity
    line_start_velocity = 0;
    while numel(regexp(line,pattern_velocity)) < 1
        line_start_velocity = line_start_velocity + 1;
        line = fgetl(fid);
    end

    frewind(fid)
    VEL = textscan(fid,'%f%f%f','HeaderLines',line_start_velocity);
    VEL = cell2mat(VEL);
        
    fclose(fid);
    
    UX = reshape(VEL(:,1),Nx,Ny,Nz);
    UY = reshape(VEL(:,2),Nx,Ny,Nz);
    UZ = reshape(VEL(:,3),Nx,Ny,Nz);
    
    slice = 101; % Middle plane

    UX_slice = UX(:,:,slice)';
    UY_slice = UY(:,:,slice)';
    UZ_slice = UZ(:,:,slice)';

    NORM = sqrt(UX_slice.^2 + UY_slice.^2 + UZ_slice.^2);
    u_max_LBM = max(max(NORM));

end

% Parameters of the simulation
alpha = 35;
a_fib = 3e-7;
B = 8e-25;
R_cyl = 8.3e-6;
dx = 120;
dy = 120;
lambda = 30e-6;
mu = 6.1e-3;
u_max_phys = 200e-6;
scale = 0.25e-6;

% Mechanical properties of the fiber
S = 4*B/a_fib^2;
ks = 100*S/(2*a_fib);
kb = B/(2*a_fib);

% Position of the pillars
nx = 120;
ny = 120;
nz = 200;
n_obs = 4;

r_obs(1,:) = [0, 0];
r_obs(2,:) = [0, ny];
r_obs(3,:) = [nx, ny];
r_obs(4,:) = [nx, 0];

r_obs = r_obs * scale;

% Data save location
data_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Simulation\Simulations_LubricationModel'];

% Figure location
figure_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Simulation\Gamma analysis'];

% Data location
parent_path = 'E:\Experimental Data (RAW)\Simulation\data';
sub1_path = dir(parent_path);
sub1_path = sub1_path(arrayfun(@(x) startsWith(x.name, 'n'), sub1_path));
average_gamma_for_one_length = zeros(1,length(sub1_path));
average_gamma_for_one_length_signConsidered = zeros(1,length(sub1_path));
Length = zeros(1,length(sub1_path));
for sub1Path_i = 1:length(sub1_path)
    current_beadsNum = sub1_path(sub1Path_i).name;

    sub2_path = dir(fullfile(parent_path, current_beadsNum));
    sub2_path = sub2_path(arrayfun(@(x) startsWith(x.name, 'simu'), sub2_path));

    XY = cell(1,length(sub2_path)); 
    GAMMA_noTension = cell(1,length(sub2_path));
    GAMMA_SIGN_noTension = cell(1,length(sub2_path));
    average_gamma_for_one_case_noTension = zeros(1,length(sub2_path));
    average_gamma_for_one_case_signConsidered_noTension = zeros(1,length(sub2_path));
    GAMMA = cell(1,length(sub2_path));
    GAMMA_SIGN = cell(1,length(sub2_path));
    average_gamma_for_one_case = zeros(1,length(sub2_path));
    average_gamma_for_one_case_signConsidered = zeros(1,length(sub2_path));
    for sub2Path_i = 1:length(sub2_path)

        folder = fullfile(parent_path, current_beadsNum, sub2_path(sub2Path_i).name);
        num_files = length(dir(fullfile(folder, 'output_data', '*.vtk')));
        disp(['Number of files in output_data: ', num2str(num_files)]);

        vec_n = 2:num_files; % The iterations you want to plot (same numbers as in the file names)

        n_snaps = length(vec_n);
        colors = jet(n_snaps);

        file = strcat(folder,'/output_data/fibers_000000.vtk');
        DATA_FIBER = read_fiber(file);
        n_beads = size(DATA_FIBER,1);

        gamma = zeros(n_snaps,n_beads);
        gamma_sign = zeros(n_snaps,n_beads);

        gamma_noTension = zeros(n_snaps,n_beads);
        gamma_sign_noTension = zeros(n_snaps,n_beads);

        xy = zeros(n_snaps,n_beads);

        for i_snap = 1:n_snaps

            n = vec_n(i_snap)-1;

            % Load the beads position at the current timestep
            file = strcat(folder,'/output_data/fibers_',num2str(n,'%06i'),'.vtk');
            DATA_FIBER = read_fiber(file);

            r = [DATA_FIBER(:,1) DATA_FIBER(:,2) DATA_FIBER(:,3)];

            % Compute all the forces
            f_stretching = forces_stretching(n_beads,r,ks,a_fib);
            f_bending = forces_bending(n_beads,r,kb);
            f_steric_fib_fib = forces_steric_fiber_fiber(n_beads,a_fib,r,mu,u_max_phys,1.1,1);
            f_steric_fib_obs = forces_steric_fiber_obstacles_cylinders(nx,ny,nz,n_beads,n_obs,a_fib,r,r_obs,R_cyl,mu,u_max_phys,1.1,1,scale);

            F = f_stretching + f_bending + f_steric_fib_fib + f_steric_fib_obs;
            F_noTension = f_bending + f_steric_fib_fib + f_steric_fib_obs;

            % Compute U = MF
            radi_vec = a_fib*ones(n_beads,1);
            M = compute_mobility_matrix(n_beads,r,radi_vec,mu);
            u_fib = reshape_vectors(M * reshape_inline(F));
            u_fib_noTension = reshape_vectors(M * reshape_inline(F_noTension));

            % Interpolate the LBM flow at the beads center of mass
            u_flow = lbm_interpolation_velocity(nx,ny,nz,VEL,n_beads,r,scale) / u_max_LBM*u_max_phys;

            % Add the contribution of the flow to the fiber beads velocity
            u_fib = u_fib + u_flow;
            u_fib_noTension = u_fib_noTension + u_flow;

            % Add the contribution of the lubrication force
            f_lubrication = forces_lubrication(nx,ny,nz,n_beads,n_obs,a_fib,r,r_obs,R_cyl,mu,u_fib,scale);
            f_lubrication_noTension = forces_lubrication(nx,ny,nz,n_beads,n_obs,a_fib,r,r_obs,R_cyl,mu,u_fib_noTension,scale);

            u_fib = u_fib + reshape_vectors(M * reshape_inline(f_lubrication));
            u_fib_noTension = u_fib_noTension + reshape_vectors(M * reshape_inline(f_lubrication_noTension));

            if boolPlotVelocities

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % PLOT THE VELOCITY OF THE FIBER BEADS AS VECTORS
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                scale_vectors = 5e-2;
                angle_arrows = 30;
                scale_arrows = 0.3;

                figure, hold on;

                plot(r(:,1),r(:,2),'Color',colors(i_snap,:),'LineWidth',2);

                for i = 1:4:n_beads

                    color = 'r';

                    x0 = r(i,1);
                    y0 = r(i,2);
                    x1 = x0 + u_fib(i,1) * scale_vectors;
                    y1 = y0 + u_fib(i,2) * scale_vectors;

                    % Plot the main segment of the vector
                    plot([x0 x1],[y0 y1],'Color',color,'LineWidth',1.5)

                    beta = atan2(y1-y0,x1-x0);
                    arrow_length = scale_arrows*sqrt((x1-x0)^2+(y1-y0)^2);

                    % First segment of the arrow
                    x2 = x1 + arrow_length*cos(beta + pi - angle_arrows*pi/180);
                    y2 = y1 + arrow_length*sin(beta + pi - angle_arrows*pi/180);
                    plot([x1 x2],[y1 y2],'Color',color,'LineWidth',1)

                    % Second segment of the arrow
                    x3 = x1 + arrow_length*cos(beta + pi + angle_arrows*pi/180);
                    y3 = y1 + arrow_length*sin(beta + pi + angle_arrows*pi/180);
                    plot([x1 x3],[y1 y3],'Color',color,'LineWidth',1)

                end


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % PLOT THE VELOCITY OF THE FLOW AS VECTORS
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                scale_vectors = 5e-2;
                angle_arrows = 30;
                scale_arrows = 0.3;

                for i = 1:4:n_beads

                    color = 'b';

                    x0 = r(i,1);
                    y0 = r(i,2);
                    x1 = x0 + u_flow(i,1) * scale_vectors;
                    y1 = y0 + u_flow(i,2) * scale_vectors;

                    % Plot the main segment of the vector
                    plot([x0 x1],[y0 y1],'Color',color,'LineWidth',1.5)

                    beta = atan2(y1-y0,x1-x0);
                    arrow_length = scale_arrows*sqrt((x1-x0)^2+(y1-y0)^2);

                    % First segment of the arrow
                    x2 = x1 + arrow_length*cos(beta + pi - angle_arrows*pi/180);
                    y2 = y1 + arrow_length*sin(beta + pi - angle_arrows*pi/180);
                    plot([x1 x2],[y1 y2],'Color',color,'LineWidth',1)

                    % Second segment of the arrow
                    x3 = x1 + arrow_length*cos(beta + pi + angle_arrows*pi/180);
                    y3 = y1 + arrow_length*sin(beta + pi + angle_arrows*pi/180);
                    plot([x1 x3],[y1 y3],'Color',color,'LineWidth',1)

                end

                axis equal

                camroll(alpha);

                set(gca,'XTick',[],'YTick',[])

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % COMPUTE THE ANGLE GAMMA
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for i = 1:n_beads
                gamma(i_snap,i) = acos( dot(u_fib(i,:),u_flow(i,:)) / (norm(u_fib(i,:)) * norm(u_flow(i,:))) ) * 180 / pi;
                the_cross_product = cross(u_fib(i,:),u_flow(i,:));
                gamma_sign(i_snap,i) = -sign(the_cross_product(3)); % positive means: u_fib counterclockwise rotates based on u_flow.

                gamma_noTension(i_snap,i) = acos( dot(u_fib_noTension(i,:),u_flow(i,:)) / (norm(u_fib_noTension(i,:)) * norm(u_flow(i,:))) ) * 180 / pi;
                the_cross_product_noTension = cross(u_fib_noTension(i,:),u_flow(i,:));
                gamma_sign_noTension(i_snap,i) = -sign(the_cross_product_noTension(3)); % positive means: u_fib_noTension counterclockwise rotates based on u_flow.
            end

            % if r(1,1) > r(end,1)
            %     gamma(i_snap,:) = flip(gamma(i_snap,:));
            %     detachment_bead = find(flip(f_steric_fib_obs(:,1)) ~= 0,1,'last');
            % else
            %     detachment_bead = find(f_steric_fib_obs(:,1) ~= 0,1,'last');
            % end
            % gamma(i_snap,1:detachment_bead) = NaN;

            attachment_bead = find(abs(f_steric_fib_obs(:,1)) > 0);
            gamma(i_snap,attachment_bead) = NaN;
            gamma_noTension(i_snap,attachment_bead) = NaN;

            xy(i_snap, :) = complex(r(:,1),r(:,2))';
            % NOTE: The transpose of a complex number makes the conjugate of the complex number, so the overall transport has been flipped!!! (Z. LI)

        end

        if boolPlotGammaProfiles
            % Plot the gamma profiles

            figure, set_plot(gcf,gca)

            for i_snap = 1:n_snaps
                plot((1:n_beads)*2*a_fib/lambda,gamma(i_snap,:),'Color',colors(i_snap,:),'LineWidth',2,'Tag',strcat('curve_n',num2str(vec_n(i_snap)))); hold on;
            end

            xlabel('$s/\lambda$')
            ylabel('$\gamma~(^\circ)$')
            xlim([0 2*n_beads*a_fib/lambda])
            ylim([0 60])
            set(gca,'XTick',0:0.4:2*n_beads*a_fib/lambda,'YTick',0:15:60)
        end

        XY{sub2Path_i} = xy;
        GAMMA{sub2Path_i} = gamma;
        GAMMA_SIGN{sub2Path_i} = gamma_sign;
        average_gamma_for_one_case(sub2Path_i) = mean(mean(gamma, 'omitnan'), 'omitnan');
        average_gamma_for_one_case_signConsidered(sub2Path_i) = mean(mean(gamma.*gamma_sign, 'omitnan'), 'omitnan');

        GAMMA_noTension{sub2Path_i} = gamma_noTension;
        GAMMA_SIGN_noTension{sub2Path_i} = gamma_sign_noTension;
        average_gamma_for_one_case_noTension(sub2Path_i) = mean(mean(gamma_noTension, 'omitnan'), 'omitnan');
        average_gamma_for_one_case_signConsidered_noTension(sub2Path_i) = mean(mean(gamma_noTension.*gamma_sign_noTension, 'omitnan'), 'omitnan');

    end

    average_gamma_for_one_length = mean(average_gamma_for_one_case, 'omitnan');
    average_gamma_for_one_length_signConsidered = mean(average_gamma_for_one_case_signConsidered, 'omitnan');

    average_gamma_for_one_length_noTension = mean(average_gamma_for_one_case_noTension, 'omitnan');
    average_gamma_for_one_length_signConsidered_noTension = mean(average_gamma_for_one_case_signConsidered_noTension, 'omitnan');

    Length = 0.6 * n_beads * 1e-6;
    % L_0(sub1Path_i) = Length;

    save([data_mother_save_path, filesep, 'FiberLength=', num2str(Length*1e6), 'um_GammaStudy.mat'], ...
        'XY', 'GAMMA', 'GAMMA_SIGN', 'GAMMA_noTension', 'GAMMA_SIGN_noTension', 'Length', ...
        'average_gamma_for_one_case', 'average_gamma_for_one_case_noTension', ...
        'average_gamma_for_one_length', 'average_gamma_for_one_length_noTension', ...
        'average_gamma_for_one_case_signConsidered', 'average_gamma_for_one_case_signConsidered_noTension',...
        'average_gamma_for_one_length_signConsidered', 'average_gamma_for_one_length_signConsidered_noTension');
end



%% Load saved gamma data and plot the averaged gamma vs. fiber lengths (possible Fig.5 in the paper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

% Geometrical properties of the array
C2C_pillar = 30 * 1e-6; % m (delta_x & delta_y)

data_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Simulation\Simulations_LubricationModel'];
saved_data = dir(data_mother_save_path);
% only keep the .mat files inluding the '_TensionStudy'
saved_data = saved_data(arrayfun(@(x) contains(x.name, '_GammaStudy.mat'), saved_data)); 

figure_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Simulation\Gamma analysis'];

L_0_group = zeros(1, length(saved_data));
average_gamma_for_one_length_group = zeros(1, length(saved_data));
average_gamma_for_one_length_signConsidered_group = zeros(1, length(saved_data));
for ii = 1:length(saved_data)

    load(fullfile(saved_data(ii).folder, saved_data(ii).name));

    L_0_group(ii) = Length;
    average_gamma_for_one_length_signConsidered_group(ii) = average_gamma_for_one_length_signConsidered;
    average_gamma_for_one_length_group(ii) = average_gamma_for_one_length;

    % clearvars XY GAMMA GAMMA_SIGN average_gamma_for_one_case average_gamma_for_one_case_signConsidered
    % close all

end

% Plot the ensemble average of gamma vs. fiber length
% Plot this after the tension-length with code: Vic_ActinPAs_Simulation_tensionStudy.m

figure('Position', [100, 100, 800, 600], 'Color', 'w');
hold on;
yyaxis right
ax = gca;  % Get current axis
ax.YColor = [166, 216, 84]/256;  % Change left y-axis color to dark green
plot(L_0_group / C2C_pillar, average_gamma_for_one_length_signConsidered_group, 'd', ...
     'MarkerFaceColor', [166, 216, 84]/256, 'MarkerEdgeColor', 'k', 'MarkerSize', 15, 'LineWidth', 1.5);
xlabel('$L/\lambda$', 'FontSize', 28, 'Interpreter', 'latex');
ylabel('$\overline{\langle\gamma{_{\rm{No\, tension}}}\rangle}$', 'FontSize', 28, 'Interpreter', 'latex');
xlim([0 3.2]); ylim([-5 10]); % grid on
set(gca, 'FontSize', 28, 'FontName', 'Times New Roman');

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',[figure_mother_save_path, filesep, ...
    'Max. tension (only contact, sorted) and ave. gamma (non-contact beads, with sign) vs. Fiber length.eps']);

saveas(gcf, [figure_mother_save_path, filesep, 'Max. tension (only contact, sorted) and ave. gamma (non-contact beads, with sign) vs. Fiber length.png']);






%% Load saved gamma data and plot the averaged gamma ON THE LAST BEAD vs. fiber lengths (possible Fig.5 in the paper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

% Geometrical properties of the array
C2C_pillar = 30 * 1e-6; % m (delta_x & delta_y)

theta = 35; % angle in degrees (counterclockwise)
Rotate_matrix = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; % rotation matrix

data_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Simulation\Simulations_LubricationModel'];
saved_data = dir(data_mother_save_path);
% only keep the .mat files inluding the '_TensionStudy'
saved_data = saved_data(arrayfun(@(x) contains(x.name, '_GammaStudy.mat'), saved_data)); 

figure_mother_save_path = ['E:\Processing & Results\Actin Filaments in Porous Media\' ...
    'Figures\Simulation\Gamma analysis'];

L_0_group = zeros(1, length(saved_data));
average_gamma_firstBead = zeros(1, length(saved_data));
average_gamma_lastBead = zeros(1, length(saved_data));
average_gamma_firstBeadafterPillar = zeros(1, length(saved_data));
for ii = 1:length(saved_data)

    load(fullfile(saved_data(ii).folder, saved_data(ii).name));

    XY = cellfun(@(x) x', XY, 'UniformOutput', false);
    % NOTE: To transpose back!!! (Z. LI)

    L_0_group(ii) = Length;

    gamma_firstBead = [];
    gamma_lastBead = [];
    gamma_firstBeadafterPillar = [];
    for case_no = 1: length(GAMMA)
        XY_first = [real(XY{1, case_no}(1, :)); imag(XY{1, case_no}(1, :))];
        XY_last = [real(XY{1, case_no}(end, :)); imag(XY{1, case_no}(end, :))];

        XY_first_rotated = Rotate_matrix * XY_first;
        XY_last_rotated = Rotate_matrix * XY_last;

        If_inOrder = XY_last_rotated(1, :) >= XY_first_rotated(1, :);

        gamma = GAMMA{1, case_no}';
        gamma_sign = GAMMA_SIGN{1, case_no}';
        gamma = gamma_sign .* gamma;
        gamma(:, ~If_inOrder) = flipud(gamma(:, ~If_inOrder));
        gamma(:, ~any(isnan(gamma), 1)) = []; % Only keep the contact cases

        % Find the first non-NaN value after NaN for each column
        for col = 1:size(gamma, 2)
            nan_indices = find(isnan(gamma(:, col)));         
            first_non_nan_index = find(~isnan(gamma(nan_indices(end)+1:end, col)), 1, 'first') + nan_indices(end);
            if ~isempty(first_non_nan_index)
                gamma_firstBeadafterPillar = [gamma_firstBeadafterPillar, gamma(first_non_nan_index, col)];
            else 
                gamma_firstBeadafterPillar = [gamma_firstBeadafterPillar, NaN];                
            end
        end

        gamma_firstBead = [gamma_firstBead, gamma(1, :)];
        gamma_lastBead = [gamma_lastBead, gamma(end, :)];

        % % % % Plot the selected fiber
        % % % figure("Color", "w");
        % % % X_pillar = 0:C2C_pillar:1e-4; Y_pillar = 0:-C2C_pillar:-1e-4;
        % % % [X_pillar, Y_pillar] = meshgrid(X_pillar, Y_pillar);
        % % % viscircles([X_pillar(:), Y_pillar(:)], 8.3*1e-6*ones(length(X_pillar(:)), 1), ...
        % % %     'Color', [0.3 0.3 0.3], 'LineWidth', 0.3, 'EnhanceVisibility', false);
        % % % hold on
        % % % % Pick one fiber to plot
        % % % for fooo = 1:20
        % % %     x_coords = real(XY{1, case_no}(:, fooo));
        % % %     y_coords = imag(XY{1, case_no}(:, fooo));
        % % %     if If_inOrder(fooo)
        % % %         plot(x_coords, y_coords, '-om', 'LineWidth', 0.5, 'MarkerSize', 2);
        % % %     else
        % % %         x_coords = flipud(x_coords); y_coords = flipud(y_coords);
        % % %         plot(x_coords, y_coords, '-om', 'LineWidth', 0.5, 'MarkerSize', 2);
        % % %     end
        % % %     % Highlight the last bead with a different symbol
        % % %     hold on;
        % % %     plot(x_coords(end), y_coords(end), 's', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
        % % %     axis equal;
        % % % end
    end

    average_gamma_firstBead(ii) = mean(gamma_firstBead, 'omitnan');
    average_gamma_lastBead(ii) = mean(gamma_lastBead, 'omitnan');
    average_gamma_firstBeadafterPillar(ii) = mean(gamma_firstBeadafterPillar, 'omitnan');

end

% Plot the ensemble average of gamma ON THE LAST BEAD vs. fiber length
% Plot this after the tension-length with code: Vic_ActinPAs_Simulation_tensionStudy.m

myColor = cmocean('balance', 10); % colormap

figure('Position', [100, 100, 1200, 600], 'Color', 'w');
% hold on;
% yyaxis right
% ax = gca;  % Get current axis
% ax.YColor = [126, 156, 64]/256;  % Change left y-axis color to dark green
plot(L_0_group / C2C_pillar, average_gamma_firstBead, 'o', ...
     'MarkerFaceColor', myColor(3, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 15, 'LineWidth', 1.5);
hold on;
plot(L_0_group / C2C_pillar, average_gamma_lastBead, 'o', ...
     'MarkerFaceColor', myColor(6, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 15, 'LineWidth', 1.5);
hold on;
plot(L_0_group / C2C_pillar, average_gamma_firstBeadafterPillar, 'o', ...
     'MarkerFaceColor', myColor(9, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 15, 'LineWidth', 1.5);
xlabel('$L/\lambda$', 'FontSize', 28, 'Interpreter', 'latex');
ylabel('$\langle\gamma\rangle$', 'FontSize', 28, 'Interpreter', 'latex');
legend({'${\rm{First\, bead}}$', '${\rm{Last\, bead}}$','${\rm{First\,dsetached\,bead}}$'}, ...
    'Interpreter', 'latex', 'FontSize', 24, 'Location', 'bestoutside', 'Box', 'off');
xlim([0 3.2]); ylim([-40 80]); grid on
set(gca, 'FontSize', 28, 'FontName', 'Times New Roman');

set(gcf,'renderer','Painters');
print('-depsc2','-tiff','-r100','-vector',[figure_mother_save_path, filesep, ...
    'ave. gamma (last bead, with sign) vs. Fiber length.eps']);

saveas(gcf, [figure_mother_save_path, filesep, 'ave. gamma vs. Fiber length.png']);















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the position of the fiber beads

function X = read_fiber(file)
    pattern_start_points = 'POINTS';
    pattern_scientific_notation = '(+|-)?\d+(\.\d+)?(E(+|-)?\d+)?';
    fid = fopen(file,'r');

    line = fgetl(fid);

    while numel(regexp(line,pattern_start_points)) < 1
        line = fgetl(fid);
    end

    n_points = str2double(regexp(line,'\d*','match'));

    X = zeros(n_points,3);

    for i = 1:n_points
        line = fgetl(fid);
        numbers = str2double(regexp(line,pattern_scientific_notation,'match'));
        X(i,1) = numbers(1);
        X(i,2) = numbers(2);
        X(i,3) = numbers(3);        
    end

    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the stretching force

function fs = forces_stretching(n_beads,r,ks,a)

    fs = zeros(n_beads,3);

    for i = 1:n_beads

        R_i_ip_norm = 0;
        R_i_im_norm = 0;

        if (i < n_beads)
            R_i_ip_vec = r(i,:) - r(i+1,:);
            R_i_ip_norm = norm(R_i_ip_vec);
            fs_i_ip = -ks*(R_i_ip_norm - 2*a)*R_i_ip_vec/R_i_ip_norm;
        else
            fs_i_ip = 0;
        end

        if (i > 1)
            R_i_im_vec = r(i-1,:) - r(i,:);
            R_i_im_norm = norm(R_i_im_vec);

            fs_i_im = -ks*(R_i_im_norm - 2*a)*(-R_i_im_vec/R_i_im_norm);
        else
            fs_i_im = 0;
        end

        fs(i,:) = fs_i_ip + fs_i_im;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the bending force

function fb = forces_bending(n_beads,r,kb)

    fb = zeros(n_beads,3);

    % First point
    t2 = r(1,:) - r(2,:);
    t3 = r(2,:) - r(3,:);

    l2 = norm(t2);
    l3 = norm(t3);

    fb(1,:) = kb * (t3/l2/l3 + dot(t2,t3)*(-t2)/l2^3/l3);

    % Second point
    t2 = r(1,:) - r(2,:);
    t3 = r(2,:) - r(3,:);
    t4 = r(3,:) - r(4,:);

    l2 = norm(t2);
    l3 = norm(t3);
    l4 = norm(t4);

    term1 = (t2 - t3)/l2/l3 + dot(t2,t3)*(t2/l2^3/l3 - t3/l2/l3^3);
    term2 = t4/l3/l4 + dot(t3,t4)*(-t3)/l3^3.0d0/l4;

    fb(2,:) = kb * (term1 + term2);

    % Points 3 to n_beads-2
    for i = 3:n_beads-2

        tim  = r(i-2,:) - r(i-1,:);
        ti   = r(i-1,:) - r(i,:);
        tip  = r(i,:)   - r(i+1,:);
        tip2 = r(i+1,:) - r(i+2,:);

        lim  = norm(tim);
        li   = norm(ti);
        lip  = norm(tip);
        lip2 = norm(tip2);

        term1 = -tim/lim/li + dot(tim,ti)*ti/lim/li^3;
        term2 = (ti - tip)/li/lip + dot(ti,tip)*(ti/li^3/lip - tip/li/lip^3);
        term3 = tip2/lip/lip2 + dot(tip,tip2)*(-tip)/lip^3/lip2;

        fb(i,:) = kb * (term1 + term2 + term3);

    end

    % Point n_beads-1
    tim = r(n_beads-3,:) - r(n_beads-2,:);
    ti  = r(n_beads-2,:) - r(n_beads-1,:);
    tip = r(n_beads-1,:) - r(n_beads,:);

    lim = norm(tim);
    li  = norm(ti);
    lip = norm(tip);

    term1 = -tim/lim/li + dot(tim,ti)*ti/lim/li^3;
    term2 = (ti - tip)/li/lip + dot(ti,tip)*(ti/li^3/lip - tip/li/lip^3);

    fb(n_beads-1,:) = kb * (term1 + term2);

    % Point n_beads
    tim = r(n_beads-1,:) - r(n_beads-2,:);
    ti =  r(n_beads,:)   - r(n_beads-1,:);

    lim = norm(tim);
    li  = norm(ti);

    tim_hat  = tim/lim;
    ti_hat   = ti/li;

    term1 = 1/li*tim_hat;
    term2 = -(1/li*dot(tim_hat,ti_hat))*ti_hat;

    fb(n_beads,:) = kb * (term1 + term2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the steric force between the fiber beads

function f_steric = forces_steric_fiber_fiber(n_beads,a_fib,r_fib,mu,u_max,khi,gamma)

    f_steric = zeros(n_beads,3);

    for i = 1:n_beads
        for j = i+1:n_beads
            Rij_vec = r_fib(i,:) - r_fib(j,:);
            Rij = norm(Rij_vec);

            R_ref = khi*(a_fib + a_fib);

            if (Rij < R_ref)

                DF = 6*pi*mu*a_fib*u_max;

                xi_0 = R_ref^2 - Rij^2;
                xi_1 = R_ref^2 - (a_fib + a_fib)^2;

                strength = DF/(a_fib + a_fib);

                f_steric(i,1) = f_steric(i,1) + strength*((xi_0/xi_1)^(2*gamma))*Rij_vec(1);
                f_steric(i,2) = f_steric(i,2) + strength*((xi_0/xi_1)^(2*gamma))*Rij_vec(2);
                f_steric(i,3) = f_steric(i,3) + strength*((xi_0/xi_1)^(2*gamma))*Rij_vec(3);

                f_steric(j,1) = f_steric(j,1) - strength*((xi_0/xi_1)^(2*gamma))*Rij_vec(1);
                f_steric(j,2) = f_steric(j,2) - strength*((xi_0/xi_1)^(2*gamma))*Rij_vec(2);
                f_steric(j,3) = f_steric(j,3) - strength*((xi_0/xi_1)^(2*gamma))*Rij_vec(3);

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the steric force between the fiber beads and the pillars

function f_steric = forces_steric_fiber_obstacles_cylinders(nx,ny,nz,n_beads,n_obs,a_fib,r_fib,r_obs,R_cyl,mu,u_max,khi,gamma,scale)

    R_ref = R_cyl + khi*a_fib;
    DF = 6*pi*mu*a_fib*u_max;
    strength = DF/a_fib;
    xi_1 = R_ref*R_ref - a_fib*a_fib;

    f_steric = zeros(n_beads,3);

    r_fib_mod(:,1) = mod(r_fib(:,1),nx*scale);
    r_fib_mod(:,2) = mod(r_fib(:,2),ny*scale);
    r_fib_mod(:,3) = mod(r_fib(:,3),nz*scale);

    for i = 1:n_obs
        for j = 1:n_beads

            x_fib = r_fib_mod(j,1);
            y_fib = r_fib_mod(j,2);
            z_fib = r_fib_mod(j,3);

            Rij_vec = [r_obs(i,1), r_obs(i,2), z_fib] - [x_fib, y_fib, z_fib];
            Rij = norm(Rij_vec);

            if (Rij < R_ref)

                xi_0 = R_ref*R_ref - Rij*Rij;
                f_steric(j,:) = f_steric(j,:) - strength*((xi_0/xi_1)^(2*gamma))*Rij_vec;

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the lubrication force

function f_lubrication = forces_lubrication(nx,ny,nz,n_beads,n_obs,a_fib,r_fib,r_obs,R_cyl,mu,u_fib,scale)

    f_lubrication = zeros(n_beads,3);

    r_fib_mod(:,1) = mod(r_fib(:,1),nx*scale);
    r_fib_mod(:,2) = mod(r_fib(:,2),ny*scale);
    r_fib_mod(:,3) = mod(r_fib(:,3),nz*scale);

    f_lub_n = zeros(1,3);
    f_lub_t = zeros(1,3);

    h_min = 0.8*a_fib;

    for i = 1:n_obs
        for j = 1:n_beads

            x_fib = r_fib_mod(j,1);
            y_fib = r_fib_mod(j,2);
            z_fib = r_fib_mod(j,3);

            Rij_vec = [r_obs(i,1), r_obs(i,2), z_fib] - [x_fib, y_fib, z_fib];
            Rij = norm(Rij_vec);

            delta = Rij - (a_fib + R_cyl);
            h = max(delta,h_min);

            unit_vec_norm = -Rij_vec/Rij;
            unit_vec_tang = [-unit_vec_norm(2), unit_vec_norm(1), 0];

            if (delta/a_fib < 1)

                f_lub_n = a_fib/h*6*pi*mu*a_fib*dot(u_fib(j,:),unit_vec_norm)*unit_vec_norm;
                f_lub_t = 16/5*pi*mu*a_fib*log(h/a_fib)*dot(u_fib(j,:),unit_vec_tang)*unit_vec_tang;

                f_lubrication(j,:) = f_lub_n + f_lub_t;

            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the mobility matrix

function M = compute_mobility_matrix(n_beads_fib,r,radi_vec,mu)

    M = zeros(3*n_beads_fib,3*n_beads_fib);
    I3 = eye(3);

    for i = 1:n_beads_fib
        for j = 1:n_beads_fib

            ai = radi_vec(i);
            aj = radi_vec(j);

            Rij_vec = r(i,:) - r(j,:);
            Rij = sum(Rij_vec.*Rij_vec)^0.5;
            outer_product_Rij = Rij_vec.*Rij_vec';

            i0 = (i - 1)*3 + 1;
            j0 = (j - 1)*3 + 1;

            if (i ~= j)
                T1 = 1 + (ai^2 + aj^2)/(3*Rij^2);
                T2 = 1 - (ai^2 + aj^2)/Rij^2;

                M(i0:i0+2, j0:j0+2) = 1/(8*pi*mu*Rij) * (T1*I3 + T2/Rij^2*outer_product_Rij);
            else
                M(i0:i0+2, j0:j0+2) = 1/(6*pi*mu*ai)*I3;
            end

        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reshape vectors inline

function v_reshaped = reshape_inline(v)

    k = 1;
    n_rows = size(v,1);
    n_cols = size(v,2);

    v_reshaped = zeros(n_rows*n_cols,1);

    for i = 1:n_rows
        for j = 1:n_cols
            v_reshaped(k) = v(i,j);
            k = k + 1;
        end 
    end 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reshape vectors multi rows

function v_reshaped = reshape_vectors(v)

    k = 1;
    n_rows = size(v,1);

    v_reshaped = zeros(n_rows/3,3);

    for i = 1:n_rows/3
        for j = 1:3
            v_reshaped(i,j) = v(k);
            k = k + 1;
        end 
    end 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert 3d indices to linear indices

function idx = convert_ijk_to_idx(nx,ny,i,j,k)

    idx = 1 + i + j*(nx+1) + k*(nx+1)*(ny+1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Interpolate the LBM velocity field at the beads COM 

function u_interp = lbm_interpolation_velocity(nx,ny,nz,u,n_beads,r_fib_in,scale)

    u_interp = zeros(n_beads,3);

    r_fib = r_fib_in / scale;

    r_fib(:,1) = mod(r_fib(:,1),nx);
    r_fib(:,2) = mod(r_fib(:,2),ny);
    r_fib(:,3) = mod(r_fib(:,3),nz);

    l = 5;

    for i_point = 1:n_beads

        u_bead = zeros(1,3);

        im2 = mod(round(r_fib(i_point,1)) - 2,nx);
        im1 = mod(round(r_fib(i_point,1)) - 1,nx);
        ip1 = mod(round(r_fib(i_point,1)) + 1,nx);
        ip2 = mod(round(r_fib(i_point,1)) + 2,nx);
        vec_i = [im2, im1, mod(round(r_fib(i_point,1)),nx), ip1, ip2];

        jm2 = mod(round(r_fib(i_point,2)) - 2,ny);
        jm1 = mod(round(r_fib(i_point,2)) - 1,ny);
        jp1 = mod(round(r_fib(i_point,2)) + 1,ny);
        jp2 = mod(round(r_fib(i_point,2)) + 2,ny);
        vec_j = [jm2, jm1, mod(round(r_fib(i_point,2)),ny), jp1, jp2];

        km2 = mod(round(r_fib(i_point,3)) - 2,nz);
        km1 = mod(round(r_fib(i_point,3)) - 1,nz);
        kp1 = mod(round(r_fib(i_point,3)) + 1,nz);
        kp2 = mod(round(r_fib(i_point,3)) + 2,nz);
        vec_k = [km2, km1, mod(round(r_fib(i_point,3)),nz), kp1, kp2];

        for idx_i = 1:l

            i = vec_i(idx_i);
            dx = abs(i - r_fib(i_point,1));

            if(dx > nx - 5)
                dx = nx - dx;
            end

            phi_x = 0;
            if (abs(dx) <= 2)
                phi_x = 0.25 * (1 + cos(0.5*pi*dx));
            end

            for idx_j = 1:l

                j = vec_j(idx_j);
                dy = abs(j - r_fib(i_point,2));

                if(dy > ny - 5)
                    dy = ny - dy;
                end

                phi_y = 0;
                if (abs(dy) <= 2)
                    phi_y = 0.25 * (1 + cos(0.5*pi*dy));
                end

                for idx_k = 1:l

                    k = vec_k(idx_k);
                    dz = abs(k - r_fib(i_point,3));
                    if(dz > nz - 5)
                        dz = nz - dz;
                    end

                    phi_z = 0;
                    if (abs(dz) <= 2)
                        phi_z = 0.25 * (1 + cos(0.5*pi*dz));
                    end

                    phi = phi_x*phi_y*phi_z;

                    ic = convert_ijk_to_idx(nx,ny,i,j,k);

                    u_bead = u_bead + phi*u(ic,:);

                end
            end
        end

        u_interp(i_point,:) = u_bead;

    end

end