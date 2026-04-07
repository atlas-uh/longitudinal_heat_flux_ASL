%%
blue = [0 0.4470 0.7410];     
orange = [0.8500 0.3250 0.0980];

fig = figure(6); clf; 
tl = tiledlayout(3,2,'TileIndexing', 'columnmajor');
tl.TileSpacing = "compact";
set(gcf,'units','normalized','OuterPosition',[0,0,0.55,0.8])

label_fontsize = 24;
gca_fontsize = 18;
legend_fontsize = 18;
m_size = 8; % marker size

nbins = 80;
k_common = logspace(log10(0.001), log10(100), nbins).';

nan_run_thresh = 1000;  
hasLongNaNRun = @(v,thresh) any(diff(find(diff([0; isnan(v); 0]))) >= thresh);
keepCol_u = ~arrayfun(@(j) hasLongNaNRun(data_sonic_SLTEST.u(:,j), nan_run_thresh), 1:size(data_sonic_SLTEST.u,2));
keepCol_v = ~arrayfun(@(j) hasLongNaNRun(data_sonic_SLTEST.v(:,j), nan_run_thresh), 1:size(data_sonic_SLTEST.v,2));
keepCol_w = ~arrayfun(@(j) hasLongNaNRun(data_sonic_SLTEST.w(:,j), nan_run_thresh), 1:size(data_sonic_SLTEST.w,2));
keepCol_T = ~arrayfun(@(j) hasLongNaNRun(data_sonic_SLTEST.T(:,j), nan_run_thresh), 1:size(data_sonic_SLTEST.T,2));

keepNeutral = abs(sonic_SLTEST.zeta) < 0.05; 
keepCol = keepCol_u & keepCol_v & keepCol_w & keepCol_T & keepNeutral; 
keepCol(388) = 0; 
keepCol(390) = 0; 

S_idx = keepCol;
zeta = sonic_SLTEST.zeta(S_idx); 

G_idx = abs(sonic_grass.zeta) < 0.05; % neutral only
ylimm = [10^-5 10^3];
xlimm = [10^-3 100];

nexttile
[fk_med, fk_q1, fk_q3] = median_curve(sonic_SLTEST.kxz(:,S_idx),sonic_SLTEST.fuu_kxz(:,S_idx),k_common);
% IQR band
ok = isfinite(fk_q1) & isfinite(fk_q3) & isfinite(k_common);
fill([k_common(ok); flipud(k_common(ok))], ...
     [fk_q3(ok); flipud(fk_q1(ok))], ...
     blue, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility','off'); hold on;
% Median line
loglog(k_common, fk_med, 'color', blue, 'LineWidth', 2);
set(gca,'fontsize',gca_fontsize,'XScale', 'log', 'XTick', logspace(-4, 2, 4),'YTick', logspace(-4, 4, 3),'XScale', 'log', 'YScale', 'log'); 

[fk_med, fk_q1, fk_q3] = median_curve(sonic_grass.kxz(:,G_idx),sonic_grass.fuu_kxz(:,G_idx),k_common);
% IQR band
ok = isfinite(fk_q1) & isfinite(fk_q3) & isfinite(k_common);
fill([k_common(ok); flipud(k_common(ok))], ...
     [fk_q3(ok); flipud(fk_q1(ok))], ...
     orange, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility','off');
% Median line
loglog(k_common, fk_med, 'color', orange, 'LineWidth', 2);

ylabel ('$ \frac{F_{uu}(k_xz)}{\sigma_u^2}$','interpreter','latex','FontWeight','bold','FontSize',label_fontsize+4,'rotation',90)
xlim(xlimm);
ylim(ylimm);
cc1 = 0.5;
kx1 = linspace(0.001,20,20); 
kxISR = linspace(0.01,500,20); 
loglog(kx1,cc1*kx1.^(-3/3),'k-','linewidth',1.5)
loglog(kxISR,cc1*kxISR.^(-5/3),'k--','linewidth',1.5)
text(0.004,200,'$-1$','fontsize',16,'interpreter','latex');
text(0.1,80,'$-5/3$','fontsize',16,'interpreter','latex')
yticks([1e-4 1e-2 1 1e2])
yticklabels({'10^{-4}','10^{-2}','10^{0}','10^{2}'})


legend('SLTEST','Grass Clearing','location','southwest','interpreter','latex')

nexttile
[fk_med, fk_q1, fk_q3] = median_curve(sonic_SLTEST.kxz(:,S_idx),sonic_SLTEST.fww_kxz(:,S_idx),k_common);
% IQR band
ok = isfinite(fk_q1) & isfinite(fk_q3) & isfinite(k_common);
fill([k_common(ok); flipud(k_common(ok))], ...
     [fk_q3(ok); flipud(fk_q1(ok))], ...
     blue, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility','off'); hold on;
% Median line
loglog(k_common, fk_med, 'color', blue, 'LineWidth', 2);
set(gca,'fontsize',gca_fontsize,'XScale', 'log', 'XTick', logspace(-4, 2, 4),'YTick', logspace(-4, 4, 3),'XScale', 'log', 'YScale', 'log'); 

[fk_med, fk_q1, fk_q3] = median_curve(sonic_grass.kxz(:,G_idx),sonic_grass.fww_kxz(:,G_idx),k_common);
% IQR band
ok = isfinite(fk_q1) & isfinite(fk_q3) & isfinite(k_common);
fill([k_common(ok); flipud(k_common(ok))], ...
     [fk_q3(ok); flipud(fk_q1(ok))], ...
     orange, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility','off');
% Median line
loglog(k_common, fk_med, 'color', orange, 'LineWidth', 2);

set(gca,'fontsize',gca_fontsize,'XScale', 'log', 'XTick', logspace(-4, 2, 4),'YTick', logspace(-4, 4, 3)); 
ylabel ('$\frac{F_{ww}(k_xz)}{\sigma_w^2}$','interpreter','latex','FontWeight','bold','FontSize',label_fontsize+4,'rotation',90)
xlim(xlimm);
ylim(ylimm);
loglog(kx1,1*kx1.^(0),'k:','linewidth',2.5)
loglog(kxISR,0.5*kxISR.^(-5/3),'k--','linewidth',1.5)
text(0.004,5,'$0$','fontsize',16,'interpreter','latex');
text(0.1,80,'$-5/3$','fontsize',16,'interpreter','latex')
yticks([1e-4 1e-2 1 1e2])
yticklabels({'10^{-4}','10^{-2}','10^{0}','10^{2}'})


nexttile
[fk_med, fk_q1, fk_q3] = median_curve(sonic_SLTEST.kxz(:,S_idx),sonic_SLTEST.fTT_kxz(:,S_idx),k_common);
% IQR band
ok = isfinite(fk_q1) & isfinite(fk_q3) & isfinite(k_common);
fill([k_common(ok); flipud(k_common(ok))], ...
     [fk_q3(ok); flipud(fk_q1(ok))], ...
     blue, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility','off'); hold on;
% Median line
loglog(k_common, fk_med, 'color', blue, 'LineWidth', 2);
set(gca,'fontsize',gca_fontsize,'XScale', 'log', 'XTick', logspace(-4, 2, 4),'YTick', logspace(-4, 4, 3),'XScale', 'log', 'YScale', 'log'); 

[fk_med, fk_q1, fk_q3] = median_curve(sonic_grass.kxz(:,G_idx),sonic_grass.fTT_kxz(:,G_idx),k_common);
% IQR band
ok = isfinite(fk_q1) & isfinite(fk_q3) & isfinite(k_common);
fill([k_common(ok); flipud(k_common(ok))], ...
     [fk_q3(ok); flipud(fk_q1(ok))], ...
     orange, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility','off');
% Median line
loglog(k_common, fk_med, 'color', orange, 'LineWidth', 2);
set(gca,'fontsize',gca_fontsize,'XScale', 'log', 'XTick', logspace(-4, 2, 4),'YTick', logspace(-4, 4, 3)); 
xlabel ('$k_x z$','interpreter','latex','FontWeight','bold','FontSize',label_fontsize-2)
ylabel ('$\frac{F_{\theta\theta}(k_xz)}{\sigma_\theta^2}$','interpreter','latex','FontWeight','bold','FontSize',label_fontsize+4,'rotation',90)
xlim(xlimm);
ylim(ylimm);
loglog(kx1,0.5*kx1.^(-3/3),'k-','linewidth',1.5)
loglog(kxISR,0.5*kxISR.^(-5/3),'k--','linewidth',1.5)
text(0.004,200,'$-1$','fontsize',16,'interpreter','latex');
text(0.1,80,'$-5/3$','fontsize',16,'interpreter','latex')
yticks([1e-4 1e-2 1 1e2])
yticklabels({'10^{-4}','10^{-2}','10^{0}','10^{2}'})



nexttile
[fk_med, fk_q1, fk_q3] = median_curve(sonic_SLTEST.kxz(:,S_idx),abs(sonic_SLTEST.fuw_kxz(:,S_idx)),k_common);
% IQR band
ok = isfinite(fk_q1) & isfinite(fk_q3) & isfinite(k_common);
fill([k_common(ok); flipud(k_common(ok))], ...
     [fk_q3(ok); flipud(fk_q1(ok))], ...
     blue, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility','off'); hold on;
% Median line
loglog(k_common, fk_med, 'color', blue, 'LineWidth', 2);
set(gca,'fontsize',gca_fontsize,'XScale', 'log', 'XTick', logspace(-4, 2, 4),'YTick', logspace(-4, 4, 3),'XScale', 'log', 'YScale', 'log'); 

[fk_med, fk_q1, fk_q3] = median_curve(sonic_grass.kxz(:,G_idx),abs(sonic_grass.fuw_kxz(:,G_idx)),k_common);
% IQR band
ok = isfinite(fk_q1) & isfinite(fk_q3) & isfinite(k_common);
fill([k_common(ok); flipud(k_common(ok))], ...
     [fk_q3(ok); flipud(fk_q1(ok))], ...
     orange, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility','off');
% Median line
loglog(k_common, fk_med, 'color', orange, 'LineWidth', 2);
set(gca,'fontsize',gca_fontsize,'XScale', 'log', 'XTick', logspace(-4, 2, 4),'YTick', logspace(-4, 4, 3)); 
ylabel ('$\left| \frac{ F_{uw}(k_xz)}{\overline{u''w''}}\right|$','interpreter','latex','Fontsize',label_fontsize+4,'rotation',90)
xlim(xlimm);
ylim(ylimm);
loglog(kx1,0.5*kx1.^(-3/3),'k-','linewidth',1.5)
loglog(kxISR,0.5*kxISR.^(-7/3),'k-.','linewidth',1.5)
text(0.004,200,'$-1$','fontsize',16,'interpreter','latex');
text(0.18,80,'$-7/3$','fontsize',16,'interpreter','latex')
yticks([1e-4 1e-2 1 1e2])
yticklabels({'10^{-4}','10^{-2}','10^{0}','10^{2}'})



nexttile
[fk_med, fk_q1, fk_q3] = median_curve(sonic_SLTEST.kxz(:,S_idx),abs(sonic_SLTEST.fwT_kxz(:,S_idx)),k_common);
% IQR band
ok = isfinite(fk_q1) & isfinite(fk_q3) & isfinite(k_common);
fill([k_common(ok); flipud(k_common(ok))], ...
     [fk_q3(ok); flipud(fk_q1(ok))], ...
     blue, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility','off'); hold on;
% Median line
loglog(k_common, fk_med, 'color', blue, 'LineWidth', 2);
set(gca,'fontsize',gca_fontsize,'XScale', 'log', 'XTick', logspace(-4, 2, 4),'YTick', logspace(-4, 4, 3),'XScale', 'log', 'YScale', 'log'); 

[fk_med, fk_q1, fk_q3] = median_curve(sonic_grass.kxz(:,G_idx),abs(sonic_grass.fwT_kxz(:,G_idx)),k_common);
% IQR band
ok = isfinite(fk_q1) & isfinite(fk_q3) & isfinite(k_common);
fill([k_common(ok); flipud(k_common(ok))], ...
     [fk_q3(ok); flipud(fk_q1(ok))], ...
     orange, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility','off');
% Median line
loglog(k_common, fk_med, 'color', orange, 'LineWidth', 2);
set(gca,'fontsize',gca_fontsize,'XScale', 'log', 'XTick', logspace(-4, 2, 4),'YTick', logspace(-4, 4, 3)); 
ylabel ('$\left| \frac{F_{w\theta}(k_xz)}{ \overline{w''\theta_v''}}\right|$','interpreter','latex','Fontsize',label_fontsize+4,'rotation',90)
xlim(xlimm);
ylim(ylimm);
loglog(kx1,0.5*kx1.^(-3/3),'k-','linewidth',1.5)
loglog(kxISR,0.5*kxISR.^(-7/3),'k-.','linewidth',1.5)
text(0.004,200,'$-1$','fontsize',16,'interpreter','latex');
text(0.18,80,'$-7/3$','fontsize',16,'interpreter','latex')
yticks([1e-4 1e-2 1 1e2])
yticklabels({'10^{-4}','10^{-2}','10^{0}','10^{2}'})



nexttile
[fk_med, fk_q1, fk_q3] = median_curve(sonic_SLTEST.kxz(:,S_idx),abs(sonic_SLTEST.fuT_kxz(:,S_idx)),k_common);
% IQR band
ok = isfinite(fk_q1) & isfinite(fk_q3) & isfinite(k_common);
fill([k_common(ok); flipud(k_common(ok))], ...
     [fk_q3(ok); flipud(fk_q1(ok))], ...
     blue, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility','off'); hold on;
% Median line
loglog(k_common, fk_med, 'color', blue, 'LineWidth', 2);
set(gca,'fontsize',gca_fontsize,'XScale', 'log', 'XTick', logspace(-4, 2, 4),'YTick', logspace(-4, 4, 3),'XScale', 'log', 'YScale', 'log'); 

[fk_med, fk_q1, fk_q3] = median_curve(sonic_grass.kxz(:,G_idx),abs(sonic_grass.fuT_kxz(:,G_idx)),k_common);
% IQR band
ok = isfinite(fk_q1) & isfinite(fk_q3) & isfinite(k_common);
fill([k_common(ok); flipud(k_common(ok))], ...
     [fk_q3(ok); flipud(fk_q1(ok))], ...
     orange, 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility','off');
% Median line
loglog(k_common, fk_med, 'color', orange, 'LineWidth', 2);
set(gca,'fontsize',gca_fontsize,'XScale', 'log', 'XTick', logspace(-4, 2, 4),'YTick', logspace(-4, 4, 3)); 
xlabel ('$k_x z$','interpreter','latex','FontWeight','bold','FontSize',label_fontsize-2)
ylabel ('$\left| \frac{F_{u\theta}(k_xz)}{\overline{u''\theta_v''}}\right|$','interpreter','latex','Fontsize',label_fontsize+4,'rotation',90)
xlim(xlimm);
ylim(ylimm);
loglog(kx1,0.5*kx1.^(-3/3),'k-','linewidth',1.5)
loglog(kxISR,0.5*kxISR.^(-7/3),'k-.','linewidth',1.5)
text(0.004,200,'$-1$','fontsize',16,'interpreter','latex');
text(0.18,80,'$-7/3$','fontsize',16,'interpreter','latex')
yticks([1e-4 1e-2 1 1e2])
yticklabels({'10^{-4}','10^{-2}','10^{0}','10^{2}'})

 exportgraphics(fig, [figures_folder 'Fig_spectra_cospectra.pdf'], 'ContentType', 'vector','BackgroundColor','none');



