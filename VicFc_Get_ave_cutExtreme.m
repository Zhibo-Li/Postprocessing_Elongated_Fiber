function L_0 = VicFc_Get_ave_cutExtreme(spl_Ls, threshold)
% the function is used to calculte the contour lenght of the filament based
% on the average without the extrame values.
%
% threshold -- threshold to truncate the extreme values

L_0_coarse = mean(spl_Ls, 'omitnan'); % roughly estimated filament length
spl_Ls(spl_Ls < L_0_coarse*(1-threshold)) = [];
spl_Ls(spl_Ls > L_0_coarse*(1+threshold)) = []; % remove the extreme value
if isempty(spl_Ls)
    L_0 = L_0_coarse;
    return
end
L_0 = mean(spl_Ls, 'omitnan');

while L_0_coarse ~= L_0 && ~isnan(L_0)
    L_0_coarse = L_0;
    spl_Ls(spl_Ls < L_0_coarse*(1-threshold)) = [];
    spl_Ls(spl_Ls > L_0_coarse*(1+threshold)) = []; % remove the extreme value
    if isempty(spl_Ls)
        L_0 = L_0_coarse;
        return
    end
    L_0 = mean(spl_Ls, 'omitnan');
end
end