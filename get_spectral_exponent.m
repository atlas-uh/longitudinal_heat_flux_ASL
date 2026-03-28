function [slope,C] = get_spectral_exponent(k,P,ind)


x = log10(k(ind));
    y = log10(abs(P(ind)));   % <-- important for cospectra

    % keep only finite values
    good = isfinite(x) & isfinite(y);
    x = x(good);
    y = y(good);

    % --- require minimum number of points
    if numel(x) < 3
        slope = NaN;
        C = NaN;
        return
    end

    % --- linear fit in log-log space
    p = polyfit(x, y, 1);

    slope = p(1);
    intercept = p(2);
    C = 10^intercept;
            


