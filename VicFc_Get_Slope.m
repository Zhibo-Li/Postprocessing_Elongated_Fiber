function theSlope = VicFc_Get_Slope(theXXXYYY, segsize, segslopethres, calclength, rmvstart, rmvend)
%% Function for slope calculation on zigzag trajectory (CoM trajectories in different flow angles).
% 
% trajectories in different flow angles 
% Input: theXXXYYY(2, n) -- the x,y coordinates of the curve.
%        segsize -- define the segment size to calculate the pre-slope.
%        segslopethres -- define the threshold to select the pre-slope.
%        calclength -- define the segment length above which could be used
%        to calculate the slope
%        rmvstart -- define the segment length to be removed on the start.
%        rmvend -- define the segment length to be removed on the end.
% Output: theSlope -- the slope.

figure;plot(theXXXYYY(1,:), theXXXYYY(2,:)); hold on

for movind = 1:size(theXXXYYY, 2)-segsize
    tmp_p = polyfit(theXXXYYY(1, movind:movind+segsize), theXXXYYY(2, movind:movind+segsize),1);
    sec_p(movind) = tmp_p(1);
end

posi_sec_p_ind = find(sec_p>segslopethres);
foo = diff(posi_sec_p_ind);
posi_sec_p_cut_ind = find(foo>1);
%         figure('Position',[-1800 0 800 600]); plot(XXX(1:end-5), sec_p);
%         figure('Position',[-900 0 800 600]); plot(XXX, YYY); axis equal; hold on
if ~isempty(posi_sec_p_cut_ind)
    for piece = 1:size(posi_sec_p_cut_ind,2)+1
        if piece == 1
            XXYY = theXXXYYY(:, posi_sec_p_ind(1):posi_sec_p_ind(posi_sec_p_cut_ind(piece)));
            if XXYY(1, end) - XXYY(1, 1) > calclength 
                [p,S] = polyfit(XXYY(1, 2:end-1), XXYY(2, 2:end-1),1);
                [y_fit,delta] = polyval(p,XXYY(1, 2:end-1),S);
                if mean(delta) < 10
                    Slope_tmp(piece) = p(1);
                else
                    Slope_tmp(piece) = nan;
                end
                plot(XXYY(1, 2:end-1), XXYY(2, 2:end-1),'r'); hold on
            end
%             XXYYY{piece} = XXYY(:, 2:end);
%             XXXYYY = XXYYY{piece};
%             plot(XXXYYY(1, :), XXXYYY(2, :),'r'); hold on
        elseif piece == size(posi_sec_p_cut_ind,2)+1
            XXYY = theXXXYYY(:, posi_sec_p_ind(posi_sec_p_cut_ind(piece-1)+1):posi_sec_p_ind(end));
            if XXYY(1, end) - XXYY(1, 1) > calclength
                [p,S] = polyfit(XXYY(1, 2:end-1), XXYY(2, 2:end-1),1);
                [y_fit,delta] = polyval(p,XXYY(1, 2:end-1),S);
                if mean(delta) < 10
                    Slope_tmp(piece) = p(1);
                else
                    Slope_tmp(piece) = nan;
                end
                plot(XXYY(1, 2:end-1), XXYY(2, 2:end-1),'r'); hold on
            end
%             XXYYY{piece} = XXYY(:, 2:end);
%             XXXYYY = XXYYY{piece};
%             plot(XXXYYY(1, :), XXXYYY(2, :),'g'); hold on
        else
            XXYY = theXXXYYY(:, posi_sec_p_ind(posi_sec_p_cut_ind(piece-1)+1):posi_sec_p_ind(posi_sec_p_cut_ind(piece)));
            if XXYY(1, end) - XXYY(1, 1) > calclength
                [p,S] = polyfit(XXYY(1, 2:end-1), XXYY(2, 2:end-1),1);
                [y_fit,delta] = polyval(p,XXYY(1, 2:end-1),S);
                if mean(delta) < 10
                    Slope_tmp(piece) = p(1);
                else
                    Slope_tmp(piece) = nan;
                end
                plot(XXYY(1, 2:end-1), XXYY(2, 2:end-1),'r'); hold on
            end
%             XXYYY{piece} = XXYY(:, 2:end);
%             XXXYYY = XXYYY{piece};
%             plot(XXXYYY(1, :), XXXYYY(2, :),'y'); hold on
        end
%         XY_plot_traj_calculatedpiece{counter} = XXXYYY; % for plotting
    end
    try
        Slope_tmp(Slope_tmp==0)=[];
        theSlope = mean(Slope_tmp,'omitnan');
    catch
        warning('No slope for this one.');
        theSlope = nan;
    end
else
    theXXXYYY(:, theXXXYYY(1, :)<rmvstart)=[]; theXXXYYY(:, theXXXYYY(1, :)>rmvend)=[];
    [p,S] =   polyfit(theXXXYYY(1, :), theXXXYYY(2, :),1);
    [y_fit,delta] = polyval(p,theXXXYYY(1, :),S);
    if mean(delta) < 10
        theSlope = p(1);
    else
        theSlope = nan;
    end
end
close
end