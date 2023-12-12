clear; close all; clc;

fiber_beads_size = 0.6; % diameter in um.
time_step = 0.01; % unit:s
B = 6.9e-26;  % Bending rigidity
PAs_C2C = 30.25; % center to center distance in X, Y direction (notice that it's not 30!)

parent_path = ['D:\Dropbox\Collaboration - LadHyX\Give_to_Zhibo_nonShared\FSI - Actin in PAs\Simulations_based_on_Experiment'];
sub_path1 = dir(parent_path);
for sub2Path_i = 3:length(sub_path1)-2
    caseName = sub_path1(sub2Path_i).name;
    newStr = strrep(caseName,'o','.');

    fileinfo = dir(fullfile(parent_path, caseName, '\output_data\*.vtk'));

    for ii = 1:length(fileinfo)
        snapshot = readVTK(fullfile(fileinfo(ii).folder, fileinfo(ii).name));
        XY = snapshot.points(:, 1:2) * 1e6; % Positions of fiber beads. (unit: um)
        xy_contact = movmean(XY, 2, 1, 'Endpoints', 'discard');
        xy_start = 2*XY(1, :) - xy_contact(1, :); xy_end = 2*XY(end, :) - xy_contact(end, :);
        xy = double([xy_start; xy_contact; xy_end]); % unit: um
        CoM_xy = mean(xy, 1); % center of mass

        spl = BSpline(xy, 'order', 2, 'nint', 10);
        seg = diff(spl);
        seg_len = sqrt(seg(:,1).^2 + seg(:,2).^2);

        [L2,R2,K2] = curvature(spl);
        R2(isnan(R2)) = inf;
        RR2 = movmean(R2,ceil(size(spl,1)/100));
        Curvature = 1./RR2; % orientation

        Energy = B / 2 * sum((Curvature(2:end)).^2 .* seg_len, 'omitnan') / 1e-6;
        L_num = sum(seg_len); % unit: um

        Gyr = 1/size(spl,1) * [sum((spl(:, 1)-CoM_xy(1)).^2),  sum((spl(:, 1)-CoM_xy(1)) .* (spl(:, 2)-CoM_xy(2)));
            sum((spl(:, 2)-CoM_xy(2)) .* (spl(:, 1)-CoM_xy(1))), sum((spl(:, 2)-CoM_xy(2)).^2)];

        [eigenV,eigenD] = eig(Gyr);
        [d,ind] = sort(diag(eigenD));
        Ds = eigenD(ind,ind);
        Vs = eigenV(:,ind);
        Lambda1 = eigenD(2,2); Lambda2 =  eigenD(1,1);

        fiber_xy_in_lattice = [mod(abs(spl(:,1)+PAs_C2C), PAs_C2C)/PAs_C2C, ...
            mod(abs(spl(:,2)-PAs_C2C), PAs_C2C)/PAs_C2C];

%         plot(spl(:,1), spl(:,2),'m.', "LineStyle","none");
%         axis equal; hold on
%         [x_plot,y_plot] = meshgrid(0:PAs_C2C:450, 0:-PAs_C2C:-450);
%         viscircles([x_plot(:), y_plot(:)], 10*ones(length(x_plot(:)), 1),'Color','r');
% 
%         plot(fiber_xy_in_lattice(:,1), fiber_xy_in_lattice(:,2),'m.', "LineStyle","none");
%         axis equal; hold on
%         viscircles([0 0; 0 1; 1 0; 1 1], 1/3*ones(4, 1),'Color','r');
%         xlim([0 1]); ylim([0 1]);
% 
%         pause(0.2)
%         close

        % fiberInfo is based on 'rotated' pillar array (the flow direction changes)
        fiberInfo(ii).name = caseName;
        fiberInfo(ii).xy = xy;
        fiberInfo(ii).spl = spl; % coordinates after B-spline interpolation
        fiberInfo(ii).length = L_num; % fiber length in um.
        fiberInfo(ii).Chi = atan(Vs(2,2)/Vs(1,2)); % orientation
        fiberInfo(ii).aniso = 1 - 4*Lambda1*Lambda2/(Lambda1+Lambda2)^2; % sphericity
        fiberInfo(ii).center = CoM_xy;
        fiberInfo(ii).Lee = norm(spl(end, :) - spl(1, :))/L_num;
        fiberInfo(ii).curvature = Curvature; % Curvature
        fiberInfo(ii).energy = Energy; % Bending energy

        fiberInfo(ii).fiber_xy_in_lattice = fiber_xy_in_lattice;
        fiberInfo(ii).fiber_CoM_xy_in_lattice = [sign(CoM_xy(:,1)) .* mod(abs(CoM_xy(:,1)), PAs_C2C)/PAs_C2C, ...
            sign(CoM_xy(:,2)) .* mod(abs(CoM_xy(:,2)), PAs_C2C)/PAs_C2C];

    end

    save_dir = [parent_path, '\simulations_results'];
    if ~exist(save_dir, "dir")
        mkdir(save_dir)
    end
    save([save_dir, filesep, caseName, '.mat'], "fiberInfo");

    clearvars fiberInfo
end



% figure('Color','white'); plot(xy(:, 1), xy(:, 2),'b*', "LineStyle","none");
% hold on
% plot(spl(:,1), spl(:,2),'r.', "LineStyle","none")
% hold on
% plot(CoM_xy(:,1), CoM_xy(:,2),'g^', "LineStyle","none")
% hold on
% scatter(spl(:,1), spl(:, 2),5, 1./RR2,"filled"); axis equal








%% FUNCTIONS
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
