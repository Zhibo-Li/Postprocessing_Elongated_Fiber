%% Persistence length from Brownian motion of actin filaments 
%  version #1 (Francesco)
%  version #2 (Zhibo)
% * read xy trajectory of filament computed thorugh elongated_objects_in_flow.m file
% * compute persistence length from the cosine correlation

%% experiment & parameters
% initialize the code
startup

% path of the experiment
basepath='D:\Dropbox\Transfer\MatlabCodes\AfterAveBGR\';

%%  experimental parameters
Obj_Mag = 0.1;  % pixel size in µm/px

% name of the mat file to read
[filename,~] = uigetfile(basepath);
load(strcat(basepath,filename),'xy','Good_case'); 

% remove the bad cases
false_ind = ~ismember(xy(1).frame, Good_case);
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

Lc = Obj_Mag.*mean(xy.arclen_spl); % mean filament length

[fin,~] = regexp(filename,'_batch');
[~,ini] = regexp(filename,'trajectory_');
savename = filename(ini+1:fin-1);

% path for the result files
pathout = strcat(basepath,'cosine-correlation\');
[status,msg,msgID] = mkdir(pathout);  
%% compute the persistence length 

N_fil = 1;  % by default
[corr_lth_bin,corr_tht_bin,corr_lth_err,corr_tht_err] = persistence_length(xy, N_fil);
corr_lth_bin = Obj_Mag*corr_lth_bin;

%% plot results

f.figure('');
scatter(corr_lth_bin,corr_tht_bin,40,[0.4667 0.6745 0.1882],'filled')
hold on
f.interp_font('latex')
f.axes('linear',[0 corr_lth_bin(end)],'log',[0 1],...
    ' $l$ $[\mu m]$ ',' $\langle cos(\theta(l)\rangle$ ',18);

[xData, yData] = prepareCurveData(corr_lth_bin, corr_tht_bin );
disp('input the data to exlude from the fit')
[exdatafrom,~] = ginput(1);
[~,ind_exclude_data_from] = min(abs(xData-exdatafrom));
ft = fittype( 'exp(-x/(2*lp))', 'independent', 'x', 'dependent', 'y' );
excludedPoints = excludedata( xData, yData, 'Indices',(ind_exclude_data_from:1:length(xData)));
ftopts = fitoptions( 'Method', 'NonlinearLeastSquares' );
ftopts.Display = 'Off';
ftopts.StartPoint = 15;
ftopts.Exclude = excludedPoints;
ftopts.Weights = 1./corr_tht_err;

close all
%% Fit model to data and replot.
[fitres, GOF, j] = fit(xData, yData, ft, ftopts) %#ok<NOPTS>
Lp = fitres.lp;
ci = confint(fitres);
deltaLp = ci(2)-mean(ci);

f.figure('');
scatter(corr_lth_bin(1:ind_exclude_data_from),....
    corr_tht_bin(1:ind_exclude_data_from),40,[0.4667 0.6745  0.1882],'filled')
hold on
X = linspace(0,corr_lth_bin(ind_exclude_data_from),400);
plot(X, feval(fitres,X),'--k','linewidth',1.5)
f.interp_font('latex')
f.axes('linear',[0 exdatafrom+5],'log',[0 1.05],...
    ' $l$ $[\mu m]$ ',' $\langle cos(\theta(l)\rangle$ ',18);

%% save results

file2save = strcat(pathout,'lp_',savename);
save(strcat(file2save,'.mat'),'corr_lth_bin', 'corr_lth_err', 'corr_tht_bin',...
    'corr_tht_err','Lc','N_fil','xy','ft','ftopts','fitres','Lp','deltaLp')

figHandles = findall(0,'Type','figure'); 
savefig(figHandles,strcat(file2save,'.fig'))




function [clb,ctb,clbe,ctbe] = persistence_length(xy, N_fil)
%% calculation of the persistence length

n = cell(N_fil,1);
for i = 1 : N_fil
    if  isempty(xy(i).spl)~=1
        n{i} = 0.5*cellfun(@numel,xy(i).spl)-1; 
        % number of segments for the B-spline centerline (i) in all frames
    end
end

lsegment = cell(N_fil,1);
for i = 1 : N_fil
    if  isempty(xy(i).spl)~=1
        lsegment{i} = cellfun(@(x) mean(x, 'omitnan'),xy(i).seglen_spl); 
        % typical length of a segment in the B-spline centreline
    end
end

%% compute the cosine correlation function

for i = 1: N_fil
    for j = 1 : xy(i).nframe

        tnt = [diff(xy(i).spl{j}(:,1)),diff(xy(i).spl{j}(:,2))]; % tangent to the centerline at each point
        tlength = 0; % set length counter to 0
        xlength = 0; % set length counter to 0
        P_len = [];  % set persistence length counter to empty

        for l_step = 1 : n{i}(j) - 1  % l_step is the length-step in the computation of the correlation function
            P_len(l_step) = 0; % set the counter
            l_length = movsum( xy(i).seglen_spl{j},[0 l_step-1]); % moving sum towards the right, @ each point, with window set by l
            tlength(l_step) = mean( l_length(1:end-l_step) ); % average over every starting point
            xlength(l_step) = l_step*lsegment{i}(j); % another way to compute the length l.
            for h = 1 : n{i}(j) - l_step
                % compute the cosine correlation sum
                costheta = (dot(tnt(h,:),tnt(h+l_step,:))/(norm(tnt(h,:))*norm(tnt(h+l_step,:))));
                P_len(l_step)  = P_len(l_step) +  costheta / (n{i}(j)-l_step);
            end
        end

        correlation_theta{i,j} = P_len;
        correlation_length{i,j} = tlength;  
        correlation_length2{i,j} = xlength;

    end
end

%% Average over different filaments in the FOV & all the frames in the movie

% combine all previous results
for i = 1 : N_fil
    c_tht{i} = cell2mat(horzcat(correlation_theta(i,:)))';
    c_lth{i} = cell2mat(horzcat(correlation_length(i,:)))';
end
corr_lth=cell2mat(vertcat(c_lth(:)));
corr_tht=cell2mat(vertcat(c_tht(:)));

% sort the results in ascending order
[corr_lth,index_crr] = sort(corr_lth);
corr_tht = corr_tht(index_crr);

% data binning
pitch = mean(cellfun(@mean, lsegment), 'omitnan'); % mean pitch of the bin
edges = min(corr_lth)-1/2*pitch : pitch : max(corr_lth)+1/2*pitch; % array of bin where to count
subs_array = discretize(corr_lth,edges); % discretize the data such as every bin in corr_lth has an identical number

% results of the binning operation
try
    clb = splitapply(@mean,corr_lth,subs_array);
    ctb = splitapply(@mean,corr_tht,subs_array);
    clbe = splitapply(@std,corr_lth,subs_array);
    ctbe = splitapply(@std,corr_tht,subs_array);

catch
    subs_array = subs_array(1: round(0.97* numel(subs_array)));
    corr_tht = corr_tht(1: round(0.97* numel(corr_tht)));
    corr_lth = corr_lth(1: round(0.97* numel(corr_lth)));

    clb = splitapply(@mean,corr_lth,subs_array);
    ctb = splitapply(@mean,corr_tht,subs_array);
    clbe = splitapply(@std,corr_lth,subs_array);
    ctbe = splitapply(@std,corr_tht,subs_array);
end

% remove nan indexes
% el = isnan(subs_array);
% subs_array(el==1)=[];
% corr_lth(el==1)=[];
% corr_tht(el==1)=[];

% clb = accumarray(subs_array,corr_lth,[],@mean);
% ctb = accumarray(subs_array,corr_tht,[],@mean);
% clbe = accumarray(subs_array,corr_lth,[],@std);
% ctbe = accumarray(subs_array,corr_tht,[],@std);

end



