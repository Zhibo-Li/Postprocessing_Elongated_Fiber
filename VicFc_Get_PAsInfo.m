function [PA_Ra, PA_cl, PA_rw] = VicFc_Get_PAsInfo(PAs)
%% Function for pillar array information.
% 
% To extract the pillar columns and rows.
% PAs: the path of the PAs results got from the 'Circle Finder' APP.
% Including variables 'centers', 'circleMask', 'metric', and 'radii'.

load(PAs);  % Load the PAs information.
PA_Ra = mean(radii);  % The averaged the radius.
sorted_centers = sort(centers(:,1));  % Calculate the interval along x-direction.
Jumps = find(diff(sorted_centers) > 100); Jumps = [0;Jumps;size(sorted_centers,1)];  % Please change the value '100' accordingly.
for ii = 1:size(Jumps, 1)-1
    PA_cl(ii) = mean(sorted_centers(Jumps(ii)+1: Jumps(ii + 1)));
end

sorted_centersY = sort(centers(:,2));  % Calculate the interval along y-direction.
JumpsY = find(diff(sorted_centersY) > 80); JumpsY = [0;JumpsY;size(sorted_centersY,1)];   % Please change the value '80' accordingly.
for ii = 1:size(JumpsY, 1)-1
    PA_rw(ii) = mean(sorted_centersY(JumpsY(ii)+1: JumpsY(ii + 1)));
end

end