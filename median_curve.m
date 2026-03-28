function [fk_med, fk_q1, fk_q3] = median_curve(k, fk, k_common)
% Interpolate each run onto k_common in log-space, then compute
% median and interquartile range (25th/75th percentiles) across runs.

nbins = length(k_common);
nRuns = size(k,2);

Smat = nan(nbins, nRuns);

for j = 1:nRuns
    kj = k(:,j);
    fj = fk(:,j);

    ok = isfinite(kj) & isfinite(fj) & (kj > 0);
    kj = kj(ok);
    fj = fj(ok);

    % Need at least 2 points to interpolate
    if numel(kj) < 2
        continue
    end

    % Ensure monotonic k for interp1
    [kj, ii] = sort(kj);
    fj = fj(ii);

    % Interpolate in log(k) to your common log-spaced grid
    Smat(:,j) = interp1(log(kj), fj, log(k_common), 'linear', NaN);
end

fk_med = median(Smat, 2, 'omitnan');
fk_q1  = prctile(Smat, 10, 2);  
fk_q3  = prctile(Smat, 90, 2);   
