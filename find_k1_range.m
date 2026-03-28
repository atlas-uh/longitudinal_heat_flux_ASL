function ind = find_k1_range(k, Pk, win, slope_tol, kz_min, kz_max, z)
    logk = log10(k(:));
    logPk = log10(abs(Pk(:)));
    
    local_slope = NaN(size(k(:)));
    hw = floor(win/2);
    for j = (hw+1):(length(k)-hw)
        idx = (j-hw):(j+hw);
        p = polyfit(logk(idx), logPk(idx), 1);
        local_slope(j) = p(1);
    end
    
    candidates = local_slope > (-1 - slope_tol) & ...
                 local_slope < (-1 + slope_tol) & ...
                 ~isnan(local_slope);
    
    % Longest contiguous block
    d = diff([0; candidates(:); 0]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    
    if isempty(starts) || max(ends - starts) < 3
        % Fallback: use manual kz range
        ind = (k.*z > kz_min) & (k.*z < kz_max);
        return
    end
    
    [~, longest] = max(ends - starts);
    ind = false(size(k(:)));
    ind(starts(longest):ends(longest)) = true;
    ind = ind(:)';
end