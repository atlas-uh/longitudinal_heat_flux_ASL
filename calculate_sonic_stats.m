function [output] = calculate_sonic_stats(sonic_struct)
    nu = 1.5e-5;
    fs = sonic_struct.fs;
    z = sonic_struct.z;
    nf_frac = 0.01;  % fraction of high-freq tail for noise estimation

    for i = 1:size(sonic_struct.u, 2)

        u = sonic_struct.u(:,i);
        v = sonic_struct.v(:,i);
        w = sonic_struct.w(:,i);
        T = sonic_struct.T(:,i);

        if sum(isnan(u)) > 5 || sum(isnan(v)) > 5 || sum(isnan(w)) > 5 || sum(isnan(T)) > 5
            continue
        end

        % ----- Apply a 5-min detrending
        tcut = 5;

        [u_de,~] = detrendHutchins(u, fs, tcut);
        un = u_de - mean(u_de);

        [v_de,~] = detrendHutchins(v, fs, tcut);
        vn = v_de - mean(v_de);

        [w_de,~] = detrendHutchins(w, fs, tcut);
        wn = w_de - mean(w_de);

        [T_de,~] = detrendHutchins(T, fs, tcut);
        Tn = T_de - mean(T_de);

        
        % ----- Basic statistics
        Ubar(i) = mean(u);
        Tbar(i) = mean(T);
        Cuu(i) = var(un);
        Cvv(i) = var(vn);
        Cww(i) = var(wn);
        CTT(i) = var(Tn);
        Cuv(i) = mean(un.*vn);
        Cuw(i) = mean(un.*wn);
        Cvw(i) = mean(vn.*wn);
        CwT(i) = mean(wn.*Tn);
        CuT(i) = mean(un.*Tn);

        g = 9.81;
        kappa = 0.4;
        TKE(i) = 0.5 * (Cuu(i) + Cvv(i) + Cww(i));
        us(i) = (Cvw(i)^2 + Cuw(i)^2)^(1/4);
        Ts(i) = -CwT(i) / us(i);
        beta(i) = g / Tbar(i);
        ws(i) = (beta(i) * CwT(i) * z)^(1/3);
        L(i) = -us(i)^3 / (kappa * g * (10*eps + CwT(i) / Tbar(i)));
        zeta(i) = z / L(i);

        % ----- Spectra
        [fuw(:,i), fuT(:,i), fwT(:,i), fuu(:,i), fww(:,i), fTT(:,i), f] = compute_spectra(un, wn, Tn, fs);

        c = (2*pi / Ubar(i));
        kx(:,i) = f * c;

        % raw wavenumber auto-spectra (unnormalized)
        fuu_raw = fuu(:,i) ./ c;
        fww_raw = fww(:,i) ./ c;
        fTT_raw = fTT(:,i) ./ c;

        % kx-normalized auto-spectra (clean)
        fuu_kx(:,i) = fuu_raw ./ Cuu(i);
        fww_kx(:,i) = fww_raw ./ Cww(i);
        fTT_kx(:,i) = fTT_raw ./ CTT(i);

        % kx-normalized cross-spectra (no noise subtraction)
        fuw_kx(:,i) = fuw(:,i) ./ c ./ Cuw(i);
        fwT_kx(:,i) = fwT(:,i) ./ c ./ CwT(i);
        fuT_kx(:,i) = fuT(:,i) ./ c ./ CuT(i);

        % kxz-normalized auto-spectra (clean)
        fuu_kxz(:,i) = fuu_raw ./ (z * Cuu(i));
        fww_kxz(:,i) = fww_raw ./ (z * Cww(i));
        fTT_kxz(:,i) = fTT_raw ./ (z * CTT(i));

        % kxz-normalized cross-spectra (no noise subtraction)
        cz = (2*pi*z / Ubar(i));
        kxz(:,i) = f * cz;
        fuw_kxz(:,i) = fuw(:,i) ./ cz ./ Cuw(i);
        fwT_kxz(:,i) = fwT(:,i) ./ cz ./ CwT(i);
        fuT_kxz(:,i) = fuT(:,i) ./ cz ./ CuT(i);

        % ----- Dissipation: isotropic (from clean unnormalized spectrum)
        eps_son.iso(i) = trapz(kx(:,i), 15*nu .* fuu_raw .* kx(:,i).^2);

        % ----- Structure functions
        clear D2 D3
        ETAs = (1/fs:1/fs:500) .* Ubar(i);
        etas = unique(floor(logspace(0, log10(length(ETAs)), 500)));
        UU = un;
        ind = 1;
        for ii = etas
            [D2(ind), D3(ind)] = structureFunctions23(UU, ii);
            ind = ind + 1;
        end
        r = ETAs(etas);
        eps_son.D3_45(i) = max(D3 / (4/5) ./ r);
        eps_son.D2(i) = max(((D2 ./ 2.12).^(3/2)) ./ r);

        lambda.D2(i) = sqrt(15*nu*Cuu(i) / eps_son.D2(i));
        Re_lambda.D2(i) = lambda.D2(i) * sqrt(Cuu(i)) / nu;

        lambda.iso(i) = sqrt(15*nu*Cuu(i) / eps_son.iso(i));
        Re = Re_lambda.D2(i);
        C3 = @(Re) (-4/5 + 8.45*(Re).^(-2/3));
        for k = 1:20
            C3(Re);
            eps_temp = max(D3 ./ C3(Re) ./ r);
            lambda_temp = sqrt(15*nu*Cuu(i) / eps_temp);
            Re = lambda_temp * sqrt(Cuu(i)) / nu;
        end
        lambda.D3(i) = lambda_temp;
        Re_lambda.D3(i) = Re;
        eps_son.D3(i) = eps_temp;
    end

    output = struct('Ubar', Ubar, 'Tbar', Tbar, ...
        'Cuu', Cuu, 'Cvv', Cvv, 'Cww', Cww, 'CTT', CTT, ...
        'Cuv', Cuv, 'Cuw', Cuw, 'CwT', CwT, 'CuT', CuT, ...
        'TKE', TKE, 'f', f, 'kx', kx, 'kxz', kxz, ...
        'fuu', fuu, 'fww', fww, 'fTT', fTT, ...
        'fuw', fuw, 'fwT', fwT, 'fuT', fuT, ...
        'fuu_kx', fuu_kx, 'fww_kx', fww_kx, 'fTT_kx', fTT_kx, ...
        'fuw_kx', fuw_kx, 'fwT_kx', fwT_kx, 'fuT_kx', fuT_kx, ...
        'fuu_kxz', fuu_kxz, 'fww_kxz', fww_kxz, 'fTT_kxz', fTT_kxz, ...
        'fuw_kxz', fuw_kxz, 'fwT_kxz', fwT_kxz, 'fuT_kxz', fuT_kxz, ...
        'beta', beta, 'us', us, 'Ts', Ts, 'ws', ws, 'z', z, 'L', L, 'zeta', zeta, ...
        'eps_son', eps_son, 'lambda', lambda, 'Re_lambda', Re_lambda);
end

function [M2, M3] = structureFunctions23(u_all, i)
    u_m = -(u_all(i+1:end,:) - u_all(1:end-i,:));
    u_m2 = (u_m).^2;
    M2 = mean(mean(u_m2));
    u_m3 = (u_m2) .* u_m;
    M3 = mean(mean(u_m3));
end