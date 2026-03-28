
function fit = fit_cospectra_CSB2(k, Pk, epsilon, A, B, gamma, ka, AuT0)
% Unified co-spectral budget solver for F_{u theta}(k)
%
% Solves the CSB ODE:
%   dF/d(lnk) + [alpha(k) + CR/AuT] * F = (1-CI)/AuT * tau_d(k) * P(k)
%
% The general solution is F = Ch*Fh + Fp, where:
%   Fh : homogeneous solution (integrating factor)
%   Fp : particular solution (variation of parameters, fully numerical)
%

k  = k(:);
Pk = Pk(:);

CR = 1.8;
CI = 3/5;

% ---- Production model ---------------------------------------------------
Prod = @(x) A ./ x .* (1 + B .* x.^2).^(-gamma);

% ---- Smooth relaxation time ---------------------------------------------
n_blend    = 4;
keff       = @(x) (x.^n_blend + ka.^n_blend).^(1/n_blend);
keff_deriv = @(x) x.^(n_blend-1) .* (x.^n_blend + ka.^n_blend).^(1/n_blend - 1);
tau_d      = @(x) epsilon.^(-1/3) .* keff(x).^(-2/3);

% ---- alpha(k) = d[ln(k/tau_d)] / d[ln k] --------------------------------
alpha_fun  = @(x) 1 + (2/3) .* x .* keff_deriv(x) ./ keff(x);

% ---- ODE grid -----------------------------------------------------------
k_min = min(k(k > 0 & isfinite(k))) * 0.5;
k_max = max(k(isfinite(k))) * 2;
N_ode = 500;
k_ode = logspace(log10(k_min), log10(k_max), N_ode)';
lnk   = log(k_ode);
dlnk  = diff(lnk);

% ---- Solve ODE for given AuT -------------------------------------------
function [Fh, Fp] = solve_ode(AuT)
    alph    = alpha_fun(k_ode);
    p_coeff = alph + CR/AuT;                          % ODE coefficient
    source  = (1-CI)/AuT .* tau_d(k_ode) .* Prod(k_ode);  % RHS

    % --- Integrating factor (shared by Fh and Fp) ---
    ln_mu = zeros(N_ode, 1);
    for i = 2:N_ode
        ln_mu(i) = ln_mu(i-1) + 0.5*(p_coeff(i-1) + p_coeff(i)) * dlnk(i-1);
    end
    mu = exp(ln_mu);

     % --- Homogeneous solution: Fh = 1/mu (unit initial condition) ---
    Fh = 1 ./ mu;

    % --- Particular solution:
    % Fp(k) = -(1/mu(k)) * integral_{k}^{k_max} mu(s)*source(s) d(lns)
    integrand = mu .* source;
    Fp = zeros(N_ode, 1);
    for i = N_ode-1:-1:1
        integrand_avg = 0.5*(integrand(i) + integrand(i+1));
        Fp(i) = (Fp(i+1)*mu(i+1) + integrand_avg*dlnk(i)) / mu(i);
    end


end

% ---- Data cleaning ------------------------------------------------------
good  = isfinite(k) & isfinite(Pk) & (k > 0) & (abs(Pk) > 0);
kfit  = k(good);
Pkfit = Pk(good);

log_Pk = log10(abs(Pkfit) + 1e-30);

% ---- For a given AuT, find optimal Ch and return residual ---------------
function [resnorm, best_c] = eval_AuT(AuT)
    [Fh_g, Fp_g] = solve_ode(AuT);

    Fh_d = interp1(lnk, Fh_g, log(kfit), 'linear', 'extrap');
    Fp_d = interp1(lnk, Fp_g, log(kfit), 'linear', 'extrap');

    % Coarse grid search for Ch
    Ch_try = logspace(-12, 12, 49);
    best_c = 1;
    best_r = Inf;
    for ic = 1:numel(Ch_try)
        F_tot = Ch_try(ic) * Fh_d + Fp_d;
        rv = sum((log10(abs(F_tot) + 1e-30) - log_Pk).^2);
        if rv < best_r
            best_r = rv;
            best_c = Ch_try(ic);
        end
    end

    % Refine Ch with fminbnd
    lo = best_c * 1e-3;
    hi = best_c * 1e3;
    try
        ref_obj = @(c) sum((log10(abs(c*Fh_d + Fp_d) + 1e-30) - log_Pk).^2);
            [best_c, best_r] = fminbnd(ref_obj, lo, hi);
    catch
    end

    resnorm = best_r;
end

% ---- Multi-start over AuT -----------------------------------------------
AuT_seeds = [AuT0, 0.5, 1, 2, 5, 10, 20, 50, 100, 500, 1000, 1e4];

best_rn_global = Inf;
best_AuT = AuT0;
best_Ch  = 1;

for j = 1:numel(AuT_seeds)
    [rn, ch] = eval_AuT(AuT_seeds(j));
    if rn < best_rn_global
        best_rn_global = rn;
        best_AuT = AuT_seeds(j);
        best_Ch  = ch;
    end
end

% Refine AuT continuously
try
    refine_AuT = @(logA) eval_AuT(exp(logA));
    [logA_opt, rn_opt] = fminbnd(refine_AuT, log(best_AuT*0.1), log(best_AuT*10));
    if rn_opt < best_rn_global
        best_AuT = exp(logA_opt);
        [~, best_Ch] = eval_AuT(best_AuT);
    end
catch
end

% ---- Final solution -----------------------------------------------------
[Fh_final, Fp_final] = solve_ode(best_AuT);

Fh_out = interp1(lnk, Fh_final, log(k), 'linear', 'extrap');
Fp_out = interp1(lnk, Fp_final, log(k), 'linear', 'extrap');

Yh   = best_Ch .* Fh_out;
Yp   = Fp_out;
Phat = Yh + Yp;

fit = struct('Ch', best_Ch, 'AuT', best_AuT, ...
             'k', k, 'Phat', Phat, 'Yh', Yh, 'Yp', Yp, ...
             'A', A, 'B', B, 'gamma', gamma, 'ka', ka);
end
