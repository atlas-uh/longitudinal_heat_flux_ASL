%% ================================================================
%  APPENDIX: Hot-wire vs Sonic validation
%  ================================================================

% --- run-index
ii_map = {[605 606], [581 582 585 586], [597 598 599]};
S = 5;  % highest hot-wire level (closest to sonic)
panelLabel = @(k) sprintf('(%c)', char('a' + (k-1)));

L_p = 0.1;  % sonic path length [m]

%% --- Figure S1: spectra and cospectra
fig_spec = figure(11); clf;
set(gcf,'units','normalized','OuterPosition',[0,0,0.6,0.75])

row_labels = {'$F_{uu}/\sigma_u^2$', '$F_{\theta\theta}/\sigma_\theta^2$', ...
              '$|F_{u\theta}/\overline{u''\theta''}|$'};

for I = 1:3
    ii = ii_map{I};
    nk = length(ii);

    % --- (row 1) velocity spectra
    subplot(3,3,I); hold on
    clear Pknorm_hw
    for n = 1:nk
        Pknorm_hw(:,n) = Up(I).Pk(:,S,n) ./ Up(I).var(:,S,n);
    end
    k_hw = Up(I).k(:,S,1);
    loglog(k_hw, mean(Pknorm_hw, 2), ...
        'color', strat_cols(S,:,I), 'linewidth', 2);
    k_son = sonic_SLTEST.kx(:,ii);
    Pk_son = sonic_SLTEST.fuu_kx(:,ii);
    loglog(mean(k_son,2), mean(Pk_son, 2), 'k-', 'linewidth', 1.5);
    set(gca, 'fontsize', 14, 'XScale', 'log', 'YScale', 'log')
    xlim([1e-3 500]); ylim([1e-4 1e2])
    if I == 1, ylabel(row_labels{1}, 'interpreter', 'latex', 'fontsize', label_fontsize); end
    title(stability_label{I}, 'fontweight', 'normal', 'fontsize', 16)
    legend('Hot-wire (1 m)', 'Sonic (2 m)', 'location', 'southwest', 'box', 'off');
    text(0.02, 0.98, panelLabel(I), 'units', 'normalized', 'fontsize', 14, 'interpreter', 'latex', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
    set(gca, 'XTick', [1e-2 1e0 1e2], 'XTickLabel', {'10^{-2}','10^{0}','10^{2}'}, 'YTick', [1e-4 1e-2 1e0 1e2],'box','on')

    xx_ls  = linspace(1e-2, 30, 20);
    xx_isr = linspace(0.1, 100, 20);
    loglog(xx_ls,  0.5*xx_ls.^(-1),   'linestyle','-', 'color',[0.7 0.7 0.7],  'linewidth', 1.5, 'HandleVisibility', 'off');
    loglog(xx_isr, 1*I*xx_isr.^(-5/3), 'linestyle','--', 'color',[0.7 0.7 0.7],'linewidth', 1.5, 'HandleVisibility', 'off');
    if I == 1
        text(0.02, 35, '$-1$', 'interpreter', 'latex', 'fontsize', 14)
        text(0.5, 10, '$-5/3$', 'interpreter', 'latex', 'fontsize', 14)
    end

    % plot k_x = 1/L_p (path-length)
    xline(1/L_p, ':', 'color', [0.5 0.5 0.5], 'linewidth', 1.2, 'HandleVisibility', 'off');
    text(1/L_p * 1.3, 0.7*min(ylim), '$1/L_p$', 'interpreter', 'latex', ...
            'fontsize', 12, 'color', 'k','HorizontalAlignment','center','VerticalAlignment','top')
  

    % --- (row 2) temperature spectra
    subplot(3,3,3+I); hold on
    clear Pknorm_hw
    for n = 1:nk
        Pknorm_hw(:,n) = Tp(I).Pk(:,S,n) ./ Tp(I).var(:,S,n);
    end
    k_hw = Tp(I).k(:,S,1);
    loglog(k_hw, mean(Pknorm_hw, 2), ...
        'color', strat_cols(S,:,I), 'linewidth', 2);
    Pk_son = sonic_SLTEST.fTT_kx(:,ii);
    loglog(mean(k_son,2), mean(Pk_son, 2), 'k-', 'linewidth', 1.5);
    set(gca, 'fontsize', 14, 'XScale', 'log', 'YScale', 'log')
    xlim([1e-3 500]); ylim([1e-4 1e2])
    if I == 1, ylabel(row_labels{2}, 'interpreter', 'latex', 'fontsize', 18); end
    text(0.02, 0.98, panelLabel(3+I), 'units', 'normalized', 'fontsize', 14, 'interpreter', 'latex', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
    set(gca, 'XTick', [1e-2 1e0 1e2], 'XTickLabel', {'10^{-2}','10^{0}','10^{2}'}, 'YTick', [1e-4 1e-2 1e0 1e2],'box','on')

    xx_ls  = linspace(1e-2, 30, 20);
    xx_isr = linspace(0.1, 100, 20);
    loglog(xx_ls,  0.5*xx_ls.^(-1),   'linestyle','-', 'color',[0.7 0.7 0.7],  'linewidth', 1.5, 'HandleVisibility', 'off');
    loglog(xx_isr, 1*I*xx_isr.^(-5/3), 'linestyle','--', 'color',[0.7 0.7 0.7],'linewidth', 1.5, 'HandleVisibility', 'off');
    if I == 1
        text(0.025, 34, '$-1$', 'interpreter', 'latex', 'fontsize', 14)
        text(0.5, 10, '$-5/3$', 'interpreter', 'latex', 'fontsize', 14)
    end

    xline(1/L_p, ':', 'color', [0.5 0.5 0.5], 'linewidth', 1.2, 'HandleVisibility', 'off');
    text(1/L_p * 1.3, 0.7*min(ylim), '$1/L_p$', 'interpreter', 'latex', ...
            'fontsize', 12, 'color', 'k','HorizontalAlignment','center','VerticalAlignment','top')
  
    % --- (row 3) uT cospectra
    subplot(3,3,6+I); hold on
    clear Pknorm_hw
    for n = 1:nk
        Pknorm_hw(:,n) = UTp(I).Pk(:,S,n) ./ UTp(I).CuT(:,S,n);
    end
    k_hw = UTp(I).k(:,S,1);
    loglog(k_hw, mean(abs(Pknorm_hw), 2), ...
        'color', strat_cols(S,:,I), 'linewidth', 2);
    Pk_son = sonic_SLTEST.fuT_kx(:,ii);
    loglog(mean(k_son,2), mean(abs(Pk_son), 2), 'k-', 'linewidth', 1.5);
    set(gca, 'fontsize', 14, 'XScale', 'log', 'YScale', 'log')
    xlim([1e-3 500]); ylim([1e-5 1e2])
    if I == 1, ylabel(row_labels{3}, 'interpreter', 'latex', 'fontsize', 18); end
    xlabel('$k_x$', 'interpreter', 'latex', 'fontsize', label_fontsize)
    text(0.02, 0.98, panelLabel(6+I), 'units', 'normalized', 'fontsize', 14, 'interpreter', 'latex', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
    set(gca, 'XTick', [1e-2 1e0 1e2], 'XTickLabel', {'10^{-2}','10^{0}','10^{2}'}, 'YTick', [1e-4 1e-2 1e0 1e2],'box','on')

    loglog(xx_ls,  0.5*xx_ls.^(-1),   'linestyle','-', 'color',[0.7 0.7 0.7],  'linewidth', 1.5, 'HandleVisibility', 'off');
    loglog(xx_isr, 2*I*xx_isr.^(-7/3), 'linestyle','-.', 'color',[0.7 0.7 0.7],'linewidth', 1.5, 'HandleVisibility', 'off');
    if I == 1
        text(0.03, 35, '$-1$', 'interpreter', 'latex', 'fontsize', 14)
        text(0.7, 10, '$-7/3$', 'interpreter', 'latex', 'fontsize', 14)
    end

    xline(1/L_p, ':', 'color', [0.5 0.5 0.5], 'linewidth', 1.2, 'HandleVisibility', 'off');
     text(1/L_p * 1.3, 0.7*min(ylim), '$1/L_p$', 'interpreter', 'latex', ...
            'fontsize', 12, 'color', 'k','HorizontalAlignment','center','VerticalAlignment','top')
  
end

exportgraphics(fig_spec, 'Fig_appendix_spectra.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');






%% --- Figure S2: sigma_u/u* and sigma_T/T* vs z/L, colored by height
fig_sT = figure; clf;
set(gcf, 'units', 'normalized', 'OuterPosition', [0, 0, 0.65, 0.35])

cmap = [158,202,225; % light blue - z1
        107,174,214; %            - z2
        66,146,198;  %            - z3
        33,113,181;  %            - z4
        8,69,148;    % dark blue  - z5
        0.00 0.00 0.00]./225;

% --- collect all data
sigu_us_all = [];
sigT_Ts_all = [];
zL_all_app  = [];
z_idx_all   = [];   % height index (1:numel(z_vec) for hot-wire, numel(z_vec)+1 for sonic)

for I = 1:3
    ii = ii_map{I};
    for n = 1:length(ii)
        us_n = sonic_SLTEST.us(ii(n));
        Ts_n = sonic_SLTEST.Ts(ii(n));
        ws_n = sonic_SLTEST.ws(ii(n));
        L_n  = sonic_SLTEST.L(ii(n));
           Tss_n = -us_n * Ts_n / ws_n; 
            uss_n = us_n^2 / ws_n;
     

        if isnan(us_n) || isnan(Ts_n) || us_n == 0 || Ts_n == 0
            continue
        end


        % hot-wire at each height
        for j = 1:numel(z_vec)
            sig_u_j = sqrt(Up(I).var(1,j,n));
            sig_T_j = sqrt(Tp(I).var(1,j,n));
            zL_j    = z_vec(j) / L_n;

            if isfinite(sig_u_j) && sig_u_j > 0
                if zL_j <0 
                    sigu_us_all = [sigu_us_all; sig_u_j / uss_n];
                    sigT_Ts_all = [sigT_Ts_all; sig_T_j / abs(Tss_n)];
                    zL_all_app  = [zL_all_app; zL_j];
                    z_idx_all   = [z_idx_all; j];
                else 
                    sigu_us_all = [sigu_us_all; sig_u_j / us_n];
                    sigT_Ts_all = [sigT_Ts_all; sig_T_j / abs(Ts_n)];
                    zL_all_app  = [zL_all_app; zL_j];
                    z_idx_all   = [z_idx_all; j];
                end 
            end
        end

        % sonic at z = 2 m
        sig_u_son = sqrt(sonic_SLTEST.Cuu(ii(n)));
        sig_T_son = sqrt(sonic_SLTEST.CTT(ii(n)));
        zL_son    = sonic_SLTEST.zeta(ii(n));

        if zL_son <0 
            sigu_us_all = [sigu_us_all; sig_u_son / uss_n];
            sigT_Ts_all = [sigT_Ts_all; sig_T_son / abs(Tss_n)];
            zL_all_app  = [zL_all_app; zL_son];
            z_idx_all   = [z_idx_all; numel(z_vec) + 1];


        else


        sigu_us_all = [sigu_us_all; sig_u_son / us_n];
        sigT_Ts_all = [sigT_Ts_all; sig_T_son / abs(Ts_n)];
        zL_all_app  = [zL_all_app; zL_son];
        z_idx_all   = [z_idx_all; numel(z_vec) + 1];
        end
    end
end

% height labels for legend
z_labels = cell(numel(z_vec) + 1, 1);
for j = 1:numel(z_vec)
    z_labels{j} = sprintf('$z=%.2f$ m', z_vec(j));
end
z_labels{end} = '$z=2$ m (sonic)';

% marker: circles for hot-wire, squares for sonic
is_son = z_idx_all == numel(z_vec) + 1;

% --- (a) sigma_u / u*
subplot(1,2,1); hold on
for j = 1:numel(z_vec)
    mask = z_idx_all == j;
    plot(-zL_all_app(mask), sigu_us_all(mask), 'o', ...
        'color', cmap(j,:), 'markerfacecolor', cmap(j,:), 'markersize', 6, ...
        'HandleVisibility', 'off');
end
plot(-zL_all_app(is_son), sigu_us_all(is_son), 's', ...
    'color', cmap(end,:), 'markerfacecolor', cmap(end,:), 'markersize', 8, ...
    'HandleVisibility', 'off');

% legend entries
for j = 1:numel(z_vec)
    plot(NaN, NaN, 'o', 'color', cmap(j,:), 'markerfacecolor', cmap(j,:), ...
        'markersize', 6, 'DisplayName', z_labels{j});
end
plot(NaN, NaN, 's', 'color', cmap(end,:), 'markerfacecolor', cmap(end,:), ...
    'markersize', 8, 'DisplayName', z_labels{end});

yline(2.5, 'k--', 'linewidth', 1.5, 'HandleVisibility', 'off');
set(gca, 'fontsize', 14,'box','on')
xlabel('$-\zeta$', 'interpreter', 'latex', 'fontsize', label_fontsize)
ylabel('$\sigma_u / u_*$', 'interpreter', 'latex', 'fontsize', label_fontsize)
legend('interpreter', 'latex', 'location', 'northwest', 'box', 'on', 'fontsize', 12)
text(0.93, 0.95, '(a)', 'units', 'normalized', 'fontsize', 14, 'interpreter', 'latex')

% --- (b) sigma_T / |T*|, colored by height, filtered
subplot(1,2,2); hold on

valid = sigT_Ts_all > 0.1;  % exclude dead T channels

for j = 1:numel(z_vec)
    mask = z_idx_all == j & valid;
    plot(-zL_all_app(mask), sigT_Ts_all(mask), 'o', ...
        'color', cmap(j,:), 'markerfacecolor', cmap(j,:), 'markersize', 6, ...
        'HandleVisibility', 'off');
end
mask_son = is_son & valid;
plot(-zL_all_app(mask_son), sigT_Ts_all(mask_son), 's', ...
    'color', cmap(end,:), 'markerfacecolor', cmap(end,:), 'markersize', 8, ...
    'HandleVisibility', 'off');

% legend entries
for j = 1:numel(z_vec)
    plot(NaN, NaN, 'o', 'color', cmap(j,:), 'markerfacecolor', cmap(j,:), ...
        'markersize', 6, 'DisplayName', z_labels{j});
end
plot(NaN, NaN, 's', 'color', cmap(end,:), 'markerfacecolor', cmap(end,:), ...
    'markersize', 8, 'DisplayName', z_labels{end});

yline(2.0, 'k--', 'linewidth', 1.5, 'HandleVisibility', 'off');
set(gca, 'fontsize', 14,'box','on')
xlabel('$-\zeta$', 'interpreter', 'latex', 'fontsize', label_fontsize)
ylabel('$\sigma_\theta / |T_*|$', 'interpreter', 'latex', 'fontsize', label_fontsize)
legend('interpreter', 'latex', 'location', 'northwest', 'box', 'on', 'fontsize', 12)
text(0.93, 0.95, '(b)', 'units', 'normalized', 'fontsize', 14, 'interpreter', 'latex')

exportgraphics(fig_sT, [figures_folder '\Fig_appendix_similarity.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');

%% --- Print near-neutral values by height
fprintf('\n=== Near-neutral (|z/L| < 0.1) ===\n')
fprintf('%-12s %10s %10s %10s %10s\n', 'Height', 'sig_u/u*', 'sig_T/T*', 'N_u', 'N_T')
fprintf('%s\n', repmat('-', 1, 55))
mask_neutral = abs(zL_all_app) < 0.1;
for j = 1:numel(z_vec)
    mj = mask_neutral & z_idx_all == j;
    mj_T = mj & valid;
    fprintf('z=%.2f m   %10.2f %10.2f %10d %10d\n', ...
        z_vec(j), mean(sigu_us_all(mj)), mean(sigT_Ts_all(mj_T)), sum(mj), sum(mj_T));
end
mj = mask_neutral & is_son;
mj_T = mj & valid;
fprintf('z=2 m (son) %9.2f %10.2f %10d %10d\n', ...
    mean(sigu_us_all(mj)), mean(sigT_Ts_all(mj_T)), sum(mj), sum(mj_T));


