%% combined production and CSB plot

% constants
CI = 3/5; 
CR = 1.8;
z_sonic = 2.0;

RGB      = orderedcolors("gem");
RGB_glow = orderedcolors("glow");

fig = figure(10); clf;
set(fig,'units','normalized','OuterPosition',[0 0 0.65 1])

tl = tiledlayout(fig, 4, 3, ...
    'TileIndexing','columnmajor', ...
    'TileSpacing','compact', ...
    'Padding','compact');

label_fontsize  = 20;
gca_fontsize    = 16;
legend_fontsize = 12;
title_fontsize  = 18;

S_list = [5 3 2];

nRows = 4;
tileIndex = @(row,col) (col-1)*nRows + row;
panelLabel = @(k) sprintf('(%c)', char('a' + (k-1)));
panelIdx = @(row,col) (row-1)*3 + col;

nf_frac = 0.1;  % top 10% for noise estimation

for I = 1:3
    if     I == 1, ii = [605 606];
    elseif I == 2, ii = [581 582 585 586];
    elseif I == 3, ii = [597 598 599];
    end

    k_son  = mean(sonic_SLTEST.kx(:,ii), 2);
    Fuw_k  = mean(sonic_SLTEST.fuw_kx(:,ii), 2);
    FwT_k  = mean(sonic_SLTEST.fwT_kx(:,ii), 2);

    for r = 1:numel(S_list)
        S = S_list(r);

        dTdz = mean(gamma(I).dTdz(S,:), 'omitnan');
        dUdz = mean(gamma(I).dUdz(S,:), 'omitnan');

        Prod_k_son_uw = Fuw_k .* dTdz;
        Prod_k_son_wT = FwT_k .* dUdz;
        Prod_k_son    = Prod_k_son_uw + Prod_k_son_wT;

        ind_son = isfinite(k_son) & isfinite(Prod_k_son);

        [k_son_uw_smooth, Prod_k_son_uw_smooth] = smooth_spectra(k_son(ind_son), Prod_k_son_uw(ind_son), 100);
        [k_son_wT_smooth, Prod_k_son_wT_smooth] = smooth_spectra(k_son(ind_son), Prod_k_son_wT(ind_son), 100);
        [k_son_smooth,    Prod_k_son_smooth]     = smooth_spectra(k_son(ind_son), Prod_k_son(ind_son),    100);

        [~,imax] = max(abs(Prod_k_son_smooth));
        ind_Prod = (k_son_smooth > k_son_smooth(imax)) & (k_son_smooth < 20) & isfinite(Prod_k_son_smooth);
        fit_Prod = fit_Production(k_son_smooth(ind_Prod), Prod_k_son_smooth(ind_Prod));


        % save for diagnostics
saved_fit_Prod(r,I).gamma = fit_Prod.gamma;
saved_fit_Prod(r,I).B     = fit_Prod.B;
saved_fit_Prod(r,I).A     = fit_Prod.A;

        % ==========================================================
        % Row 1: production plot (only for S=5)
        % ==========================================================
        if S == 5
            ax = nexttile(tl, tileIndex(1,I));
            cla(ax); hold(ax,'on');

            if ~isempty(k_son_smooth)
                loglog(ax, k_son_smooth.*z_vec(S), abs(Prod_k_son_smooth)./z_vec(S), 'k-',  'LineWidth',2,'HandleVisibility','off');
                h2 = loglog(ax, k_son_wT_smooth.*z_vec(S), abs(Prod_k_son_wT_smooth)./z_vec(S), 'k--','LineWidth',2);
                h3 = loglog(ax, k_son_uw_smooth.*z_vec(S), abs(Prod_k_son_uw_smooth)./z_vec(S), 'k:', 'LineWidth',2);
            end

            hfit = [];
            if ~isempty(fit_Prod)
                fitColor = [0 0 0];
                if exist('co','var') && size(co,1) >= 5
                    fitColor = co(5,:);
                end
                hfit = loglog(ax, fit_Prod.k.*z_vec(S), abs(fit_Prod.Pkhat)./z_vec(S), '-', ...
                    'Color', fitColor, 'LineWidth',5);
                uistack(hfit,'bottom');
            end

            set(ax,'FontSize',gca_fontsize,'YScale','log','YTick',logspace(-8,2,6));

            if I == 1
                ylabel(ax,'$|P_{u\theta}(k_xz)|$','Interpreter','latex','FontSize',label_fontsize,'Rotation',90);
                if ~isempty(hfit)
                    legend(ax,[h2,h3,hfit], ...
                        {'$\left|F_{w\theta}(k_x)\Gamma\right|$', ...
                        '$\left|F_{uw}(k_x)\Gamma_\theta\right|$', ...
                        ['$\frac{' num2str(fit_Prod.A, '%.4f') '}{k}(1+' num2str(fit_Prod.B,'%.2f') 'k^2)^{-' num2str(fit_Prod.gamma,'%.2f') '}$']}, ...
                        'Interpreter','latex','Location','southwest','FontSize',legend_fontsize+4,'Box','off');
                end
            else
                if ~isempty(hfit)
                    legend(ax,hfit, ...
                        ['$\frac{' num2str(fit_Prod.A, '%.4f') '}{k}(1+' num2str(fit_Prod.B,'%.2f') 'k^2)^{-' num2str(fit_Prod.gamma,'%.2f') '}$'], ...
                        'Interpreter','latex','Location','southwest','FontSize',legend_fontsize+4,'Box','off');
                end
            end

            title(ax,{stability_label{I};[]},'FontWeight','normal','FontSize',title_fontsize);
            text(ax,0.975,0.92, panelLabel(panelIdx(1,I)), ...
                'Units','normalized','Interpreter','latex','Color','k', ...
                'FontSize',gca_fontsize,'HorizontalAlignment','right');
            set(gca,'XScale','log','YScale','log');
            ylim([10^-5 5]);
            xlim([1e-3 100]);
            box on
        end

        % ==========================================================
        % Rows 2-4: CSB fits to hot-wire cospectra
        % ==========================================================
        row = 1 + r;

        % hot-wire cospectra (average over runs)
        k  = squeeze(mean(UTp(I).k(:,S,:), 3));
        Pk = squeeze(mean(UTp(I).Pk(:,S,:), 3));

        % smooth original for plotting (untouched)
        [k_smooth, Pk_smooth] = smooth_spectra(k, abs(Pk), 500);

        % estimate noise level for fit-range threshold
        Pk_abs = abs(Pk);
        nf_start = round(length(Pk_abs) * (1 - nf_frac));
        noise_level = median(Pk_abs(nf_start:end));

        eps_u = mean(Up(I).eps.D2(:,S,:), 'all','omitnan');
        ka = 1/z_vec(S);

        % --- Fits ---
        fit_CSB = [];
        if ~isempty(fit_Prod)
            % only fit where signal >> noise, no subtraction
            ind_fit = isfinite(k) & isfinite(Pk) & ...
                      (Pk_abs > 5*noise_level) & (k < 20/z_vec(S));
            
            [k_sm_fit, Pk_sm_fit] = smooth_spectra(k(ind_fit), Pk_abs(ind_fit), 200);

            fit_CSB = fit_cospectra_CSB2(k_sm_fit, Pk_sm_fit, eps_u, ...
                fit_Prod.A, fit_Prod.B, fit_Prod.gamma, ka, 1);

            Ch_CSB(S,I)  = fit_CSB.Ch;
            AuT_CSB(S,I) = fit_CSB.AuT;

            fprintf('I=%d S=%d (z=%.3fm): Full CSB: AuT=%.2f Ch=%.4g\n', ...
                I, S, z_vec(S), fit_CSB.AuT, fit_CSB.Ch);
        end

        % --- Plot ---
        ax = nexttile(tl, tileIndex(row,I));
        cla(ax); hold(ax,'on');

        base_color = strat_cols(4,:,I);
        light_color = min(base_color + 0.35, 1);

        % full solution
        hFull = loglog(ax, fit_CSB.k.*z_vec(S), abs(fit_CSB.Phat)./z_vec(S), '-', ...
            'LineWidth', 5, 'Color', light_color);

        % data (original, no noise subtraction)
        loglog(ax, k_smooth.*z_vec(S), abs(Pk_smooth)./z_vec(S), '-', ...
            'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');

        % components
        hHom = loglog(ax, fit_CSB.k.*z_vec(S), abs(fit_CSB.Yh)./z_vec(S), '--', ...
            'Color', base_color, 'LineWidth', 2);
        hPar = loglog(ax, fit_CSB.k.*z_vec(S), abs(fit_CSB.Yp)./z_vec(S), ':', ...
            'Color', base_color, 'LineWidth', 2);


          % AuT label
        text(ax, 0.01, 0.06, sprintf('$A_{u\\theta}=%.1f$', fit_CSB.AuT), ...
            'Units', 'normalized', 'Interpreter', 'latex', 'FontSize', legend_fontsize+4, ...
            'Color', base_color);

        xlim(ax, [1e-3 100]);
        set(ax, 'FontSize', gca_fontsize, 'XScale', 'log', 'YScale', 'log');
        box on

        if I == 1
            if S == 2
                ylabel(ax, {'$|F_{u\theta}(k_xz)|$'; '$z=0.125$ m'}, 'Interpreter','latex','FontSize',label_fontsize,'Rotation',90);
            elseif S == 3
                ylabel(ax, {'$|F_{u\theta}(k_xz)|$'; '$z=0.25$ m'},  'Interpreter','latex','FontSize',label_fontsize,'Rotation',90);
            elseif S == 5
                ylabel(ax, {'$|F_{u\theta}(k_xz)|$'; '$z=1$ m'},     'Interpreter','latex','FontSize',label_fontsize,'Rotation',90);
            end
        end

        if row == 4
            xlabel(ax, '$k_x z$', 'Interpreter','latex','FontSize',label_fontsize);
        end

        text(ax, 0.975, 0.92, panelLabel(panelIdx(row,I)), ...
            'Units','normalized','Interpreter','latex','FontSize',gca_fontsize, ...
            'HorizontalAlignment','right');

        if S == 5
            lg = legend(ax, [hFull, hHom, hPar], ...
                {'full solution', 'homogeneous', 'particular'}, ...
                'Interpreter','latex','Location','southwest', ...
                'FontSize', legend_fontsize+4, 'Box', 'off');
            lg.Position(2) = lg.Position(2) + 0.02;  % nudge up by 3% of figure height
        end

        clear hFull
    end
end

exportgraphics(fig, [figures_folder '/Fig_CSB_v2.pdf'], 'ContentType','vector', 'BackgroundColor','none');


%% --- Combined table: m and AuT by height and stability ---
stability_names = {'unstable', 'near-neutral', 'stable'};
S_list_table = [2 3 5];
z_list_table = z_vec(S_list_table);

fprintf('\n=== F_uT: ISR slope |m| and CSB parameters ===\n')
fprintf('%-14s', 'Height');
for I = 1:3, fprintf('%28s', stability_names{I}); end
fprintf('\n%-14s', '');
for I = 1:3, fprintf('  %6s %7s %7s', '|m|', 'AuT', 'Ch'); end
fprintf('\n%s\n', repmat('-', 1, 98))

for r = 1:numel(S_list_table)
    S = S_list_table(r);
    fprintf('z = %.3f m   ', z_vec(S));
    for I = 1:3
        vals = abs(m(I).phi_uT.ISR(:, S));
        vals = vals(~isnan(vals));
        if ~isempty(vals)
            m_str = sprintf('%.2f', mean(vals));
        else
            m_str = '  NaN';
        end

        if AuT_CSB(S,I) > 0
            aut_str = sprintf('%.1f', AuT_CSB(S,I));
            ch_str  = sprintf('%.3f', Ch_CSB(S,I));
        else
            aut_str = '  NaN';
            ch_str  = '  NaN';
        end

        fprintf('  %6s %7s %7s', m_str, aut_str, ch_str);
    end
    fprintf('\n');
end