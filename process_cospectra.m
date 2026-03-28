function [UTp] = process_cospectra(Up_struc, Tp_struc, pow, fs)
    S_U = find(any(isfinite(Up_struc.fluc(:,:,1)),1));
    S_T = find(any(isfinite(Tp_struc.fluc(:,:,1)),1));
    S_co = intersect(S_U, S_T);
    eps0 = 1e-12;

    for B = 1:size(Up_struc.fluc, 3)
        Tfluc_short = Tp_struc.fluc(:,S_co,B);
        Ufluc_short = Up_struc.fluc(:,S_co,B);

        UTp.CuT(:,:,B) = mean(Tfluc_short .* Ufluc_short, 1, 'omitnan');

        [temp, f] = cpsd(Ufluc_short, Tfluc_short, hamming(2^pow), 0, [], fs, 'onesided');
        UTp.Pf(:,:,B) = real(temp);
        UTp.f(:,:,B)  = f;

        Ubar = squeeze(Up_struc.mean(:,S_co,B));
        Ubar = Ubar(:).';
        UTp.Pk(:,:,B) = UTp.Pf(:,:,B) .* (Ubar/(2*pi));
        UTp.k(:,:,B)  = UTp.f(:,:,B)  .* ((2*pi)./(Ubar + eps0));

        % --- scale-wise coherence
        [Puu, ~] = pwelch(Ufluc_short, hamming(2^pow), 0, [], fs, 'onesided');
        [PTT, ~] = pwelch(Tfluc_short, hamming(2^pow), 0, [], fs, 'onesided');
        denom = sqrt(Puu) .* sqrt(PTT);
        UTp.Cs(:,:,B) = abs(UTp.Pf(:,:,B)) ./ (denom + eps0);
    end
end