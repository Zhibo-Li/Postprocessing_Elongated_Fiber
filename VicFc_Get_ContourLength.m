function [ContourL] = VicFc_Get_ContourLength(ContourL_all)
%% Function for calculating the contour length.
% 
% The contour L is the average of 10% longest snapshots except for the
% extreme values (out of 2-sigma among those averaged ones). 

Longest_10percent_num = ceil(length(ContourL_all)/10); % the number of 10% longest shapshots
if Longest_10percent_num > 3

    Sorted_L = sort(ContourL_all, 'descend');
    Longest_10percent = Sorted_L(1:Longest_10percent_num); % 10% longest L

    L_mean = mean(Longest_10percent);
    L_std = std(Longest_10percent);

    Longest_10percent(abs(Longest_10percent - L_mean) > 2*L_std) = []; % 2-sigma

    ContourL = mean(Longest_10percent);

else

    ContourL = mean(ContourL_all);

end

end