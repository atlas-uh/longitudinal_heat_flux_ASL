blue  = [0.00 0.30 0.55];
orange = [0.70 0.25 0.05];
%% ================================================================
%  FIT C1 IN -1 REGION PER RUN FOR SONIC DATA
%  ================================================================
% some restrictions
C_max = 2;      % reject |C| above this
C_min = 0.01;   % reject |C| below this
zL_max = 7;     % reject |z/L| above this

for SITE = 1:2 
    if SITE == 1 
        sonic_struc = sonic_SLTEST;
    else 
        sonic_struc = sonic_grass;
    end
    z = sonic_struc.z; 
    C_DDA_all = []; 
    C_MOST_all = []; 
    zL_all = []; 
for E = 1: length(sonic_struc.us)
    c = 2*pi / sonic_struc.Ubar(E);
    k_raw  = sonic_struc.f * c;
    Pk_raw = abs(sonic_struc.fuT(:,E) ./ c);

    % Remove bad points
    good = isfinite(k_raw) & isfinite(Pk_raw) & k_raw > 0 & Pk_raw > 0;
    if sum(good) < 20, continue; end
    k_e  = k_raw(good);
    Pk_e = Pk_raw(good);


    win = 3; 
    slope_tol = 0.75; 
    kz_min = 0.001; 
    kz_max = 1; 
    ind = find_k1_range(k_e, Pk_e, win, slope_tol, kz_min, kz_max, z);
        
    if sum(ind) < 5, continue; end

    us_e = sonic_struc.us(E);
    Ts_e = sonic_struc.Ts(E);
    L_e  = sonic_struc.L(E);
    zL_e = sonic_struc.zeta(E);

    if isnan(us_e) || isnan(Ts_e) || isnan(L_e) || us_e == 0 || Ts_e == 0
        continue
    end

    % Skip extreme stability
    if abs(zL_e) > zL_max, continue; end


    
    % MOST fit
    MOST_coeff = us_e * Ts_e;
    fit_M = fit_F1_coeff(k_e(ind), abs(Pk_e(ind)), abs(MOST_coeff));

    if abs(fit_M.C) > C_max || abs(fit_M.C) < C_min, continue; end

    C_MOST_all = [C_MOST_all; fit_M.C];
    zL_all     = [zL_all; zL_e];
   

   

    % DDA fit
    if zL_e < 0
        ws_e = real(sonic_struc.ws(E));
        if isnan(ws_e) || ws_e <= 0
            C_DDA_all = [C_DDA_all; NaN];
        else
            Tss = -us_e * Ts_e / ws_e; 
            uss = us_e^2 / ws_e;
            DDA_coeff = Tss*uss; %us_e^3 / ws_e^2 * Ts_e;
            fit_D = fit_F1_coeff(k_e(ind), abs(Pk_e(ind)), abs(DDA_coeff));
            C_DDA_all = [C_DDA_all; fit_D.C];
          end
    else
            C_DDA_all = [C_DDA_all; NaN];
    end

    % --- diagnostic plot ---
    %{
    figure(99); clf; hold on
    loglog(k_e, abs(Pk_e), '-', 'color', [0.6 0.6 0.6], 'linewidth', 1)
    loglog(k_e(ind), abs(Pk_e(ind)), '-', 'color', [0.2 0.4 0.8], 'linewidth', 2)
    k_fit = logspace(log10(min(k_e(ind))), log10(max(k_e(ind))), 50);
    loglog(k_fit, fit_M.C * abs(MOST_coeff) .* k_fit.^(-1), 'r--', 'linewidth', 1.5)
    set(gca, 'XScale', 'log', 'YScale', 'log', 'fontsize', 12, 'box', 'on')
    xlabel('$k_x$ [m$^{-1}$]', 'interpreter', 'latex', 'fontsize', 14)
    ylabel('$|F_{u\theta}|$', 'interpreter', 'latex', 'fontsize', 14)
    legend('Full spectrum', '$k^{-1}$ range', ...
        sprintf('Fit: $C_{MOST}$ = %.3f', fit_M.C), ...
        'interpreter', 'latex', 'fontsize', 12, 'box', 'off')
    title(sprintf('Site: %s  |  Run %d / %d  |  $\\zeta$ = %.3f  |  $u_* T_*$ = %.4f', ...
        sonic_labels{SITE}, E, length(sonic_struc.us), zL_e, MOST_coeff), ...
        'interpreter', 'latex', 'fontsize', 13)
    drawnow
    input('Press ENTER for next run...');
    %}

  
end


if SITE == 1 
    C_DDA_SLTEST = C_DDA_all; 
    C_MOST_SLTEST = C_MOST_all; 
    ZL_SLTEST =zL_all;
else 
    C_DDA_GRASS = C_DDA_all; 
    C_MOST_GRASS = C_MOST_all;
    ZL_GRASS =zL_all;
end


end


%%
% ================================================================
%  PLOT FIGURE: C_MOST and C_DDA vs -zeta
% ================================================================
fig7 = figure(7); clf;
set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 0.55, 0.35])

% --- (a) C_MOST ---
subplot(1,2,1); hold on

scatter(-ZL_SLTEST, C_MOST_SLTEST, 25, 'filled', 'MarkerFaceColor', blue,   'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.2, 'Marker', 'o');
scatter(-ZL_GRASS,  C_MOST_GRASS,  25, 'filled', 'MarkerFaceColor', orange, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3, 'Marker', 's');
C_unstable = mean([C_MOST_SLTEST(-ZL_SLTEST>0.5); C_MOST_GRASS(-ZL_GRASS>0.5)],'omitmissing');  
plot([0.5,zL_max],[C_unstable,C_unstable],'k:','linewidth',2); 
text(zL_max,C_unstable,[' ' num2str(C_unstable,'%.2f')],'HorizontalAlignment','left')
xlabel('$-\zeta$', 'interpreter', 'latex', 'fontsize', label_fontsize)
ylabel('$C_{MOST} = F_{u\theta}/(u_* T_* \, k_x^{-1})$', 'interpreter', 'latex', 'fontsize', label_fontsize)
xlim([0 zL_max]); ylim([0 1.5])
set(gca, 'fontsize', gca_fontsize, 'box', 'on')
text(0.01, 0.95, '(a)', 'units', 'normalized', 'interpreter', 'latex', 'fontsize', 14)

% --- inset
axes('Position',[.22 .61 .23 .27])
box on;
hold on;
scatter(-ZL_SLTEST, C_MOST_SLTEST, 25, 'filled', 'MarkerFaceColor', blue,   'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.2, 'Marker', 'o');
scatter(-ZL_GRASS,  C_MOST_GRASS,  25, 'filled', 'MarkerFaceColor', orange, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3, 'Marker', 's');

% bin-medians 
bin_edges_in = linspace(0, 0.5, 11);
bin_cen_in   = (bin_edges_in(1:end-1) + bin_edges_in(2:end)) / 2;
neg_ZL_all   = -[ZL_SLTEST; ZL_GRASS];
C_MOST_all   = [C_MOST_SLTEST; C_MOST_GRASS];
bm_MOST_in = NaN(length(bin_cen_in), 1);
for b = 1:length(bin_cen_in)
    in_b = neg_ZL_all >= bin_edges_in(b) & neg_ZL_all < bin_edges_in(b+1);
    if sum(in_b) >= 3
        errorbar(bin_cen_in(b), median(C_MOST_all(in_b)), ...
            median(C_MOST_all(in_b)) - prctile(C_MOST_all(in_b),25), ...
            prctile(C_MOST_all(in_b),75) - median(C_MOST_all(in_b)), ...
            'k-o', 'linewidth', 1, 'MarkerSize', 0.1, 'MarkerFaceColor', 'k', 'CapSize', 3, 'HandleVisibility', 'off');
        bm_MOST_in(b) = median(C_MOST_all(in_b));
    end
end
ok = ~isnan(bm_MOST_in); 
plot(bin_cen_in(ok), bm_MOST_in(ok), 'k-', 'linewidth', 1, 'HandleVisibility', 'off')
cv_bm_MOST = std(bm_MOST_in(ok)) / mean(bm_MOST_in(ok));
xlim([0 0.5]);
ylim([0 1.5]); 
set(gca,'fontsize',gca_fontsize)


% --- (b) C_DDA ---
subplot(1,2,2); hold on

scatter(-ZL_SLTEST, C_DDA_SLTEST, 25, 'filled', 'MarkerFaceColor', blue,   'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.2, 'Marker', 'o');
scatter(-ZL_GRASS,  C_DDA_GRASS,  25, 'filled', 'MarkerFaceColor', orange, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3, 'Marker', 's');
C_unstable = mean([C_DDA_SLTEST(-ZL_SLTEST>0.5&-ZL_SLTEST<zL_max); C_DDA_GRASS(-ZL_GRASS>0.5&-ZL_GRASS<zL_max)],'omitmissing');  
plot([0.5, zL_max],[C_unstable,C_unstable],'k:','linewidth',2); 
text(zL_max,C_unstable,[' ' num2str(C_unstable,'%.2f')],'HorizontalAlignment','left')
xlabel('$-\zeta$', 'interpreter', 'latex', 'fontsize', label_fontsize)
ylabel('$C_{DDA} = F_{u\theta}/(u_{**} T_{**} \, k_x^{-1})$', 'interpreter', 'latex', 'fontsize', label_fontsize)
xlim([0 zL_max]); ylim([0 1.5])
set(gca, 'fontsize', gca_fontsize, 'box', 'on')
text(0.01, 0.95, '(b)', 'units', 'normalized', 'interpreter', 'latex', 'fontsize', 14)


% --- inset
axes('Position',[.66 .61 .23 .27])
box on;
hold on;
scatter(-ZL_SLTEST, C_DDA_SLTEST, 25, 'filled', 'MarkerFaceColor', blue,   'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.2, 'Marker', 'o'); 
scatter(-ZL_GRASS,  C_DDA_GRASS,  25, 'filled', 'MarkerFaceColor', orange, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.3, 'Marker', 's');

% bin medians
C_DDA_all = [C_DDA_SLTEST; C_DDA_GRASS];
valid_DDA = ~isnan(C_DDA_all);
bm_DDA_in = NaN(length(bin_cen_in), 1);
for b = 1:length(bin_cen_in)
    in_b = neg_ZL_all >= bin_edges_in(b) & neg_ZL_all < bin_edges_in(b+1) & valid_DDA;
    if sum(in_b) >= 3
        errorbar(bin_cen_in(b), median(C_DDA_all(in_b)), ...
            median(C_DDA_all(in_b)) - prctile(C_DDA_all(in_b),25), ...
            prctile(C_DDA_all(in_b),75) - median(C_DDA_all(in_b)), ...
            'k-o', 'linewidth', 1, 'MarkerSize', 0.1, 'MarkerFaceColor', 'k', 'CapSize', 3, 'HandleVisibility', 'off');
        bm_DDA_in(b) = median(C_DDA_all(in_b));
    end
end
ok = ~isnan(bm_DDA_in);
plot(bin_cen_in(ok), bm_DDA_in(ok), 'k-', 'linewidth', 1, 'HandleVisibility', 'off')
cv_bm_DDA = std(bm_DDA_in(ok)) / mean(bm_DDA_in(ok));
xlim([0 0.5]);
ylim([0 1.5]); 
legend('SLTEST', 'Grass Clearing', 'Bin median', 'location', 'northeast', 'box', 'off', 'fontsize', 14)
set(gca,'fontsize',gca_fontsize)



exportgraphics(fig7, [figures_folder '\Fig_C_combined.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');


%% --- Print stats ---
mask_MOST_near = neg_ZL_all > 0 & neg_ZL_all < 0.5;
mask_DDA_near  = mask_MOST_near & valid_DDA;
mask_MOST_far  = neg_ZL_all > 0.5 & neg_ZL_all < zL_max;
mask_DDA_far   = mask_MOST_far & valid_DDA;

fprintf('\n=== Unstable (-zeta > 0.5) ===\n')
fprintf('C_MOST: mean=%.3f, std=%.3f, CV=%.0f%%, N=%d\n', ...
    mean(C_MOST_all(mask_MOST_far)), std(C_MOST_all(mask_MOST_far)), ...
    100*std(C_MOST_all(mask_MOST_far))/mean(C_MOST_all(mask_MOST_far)), sum(mask_MOST_far))
fprintf('C_DDA:  mean=%.3f, std=%.3f, CV=%.0f%%, N=%d\n', ...
    mean(C_DDA_all(mask_DDA_far)), std(C_DDA_all(mask_DDA_far)), ...
    100*std(C_DDA_all(mask_DDA_far))/mean(C_DDA_all(mask_DDA_far)), sum(mask_DDA_far))

fprintf('\n=== Bin-median CV (-zeta < 0.5) ===\n')
fprintf('C_MOST: mean=%.3f, std=%.3f, CV=%.0f%%\n', ...
    mean(bm_MOST_in(~isnan(bm_MOST_in))), std(bm_MOST_in(~isnan(bm_MOST_in))), ...
    100*std(bm_MOST_in(~isnan(bm_MOST_in)))/mean(bm_MOST_in(~isnan(bm_MOST_in))))
fprintf('C_DDA:  mean=%.3f, std=%.3f, CV=%.0f%%\n', ...
    mean(bm_DDA_in(~isnan(bm_DDA_in))), std(bm_DDA_in(~isnan(bm_DDA_in))), ...
    100*std(bm_DDA_in(~isnan(bm_DDA_in)))/mean(bm_DDA_in(~isnan(bm_DDA_in))))
