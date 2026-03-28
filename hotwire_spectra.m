C_MOST_all = [];
C_DDA_all  = [];
zL_all     = [];
stab_all   = [];
site_all   = {};

% coefficients for plotting purposes
coeff_AE = [0.2 0.15 0.15; ...
            0.3 0.4 0.3; ...
            0.07 0.04 0.005]; 
coeff_ISR = [0.12 0.15 0.6; ...
             0.6 0.6 0.7; ...
             0.07 0.18 0.18];

% plotting parameters 
label_fontsize = 24;
gca_fontsize = 16;
legend_fontsize = 18;

blue  = [0.00 0.30 0.55];
orange = [0.70 0.25 0.05];

%% ================================================================
%  FIGURE 1: Hot-wire spectra and cospectra (ensemble-averaged)
%  ================================================================
panelLabels = {'(a)','(b)','(c)';'(d)','(e)','(f)';'(g)','(h)','(i)'};

fig = figure(8); clf;
set(gcf,'Units','normalized','OuterPosition',[0,0,0.8,1])
for I = 1:3
    % --- ROW 1: VELOCITY SPECTRA
    subplot(3,3,I);
    S_plot = find(Up(I).Pk(1,:,1));
    for S = S_plot
        k = mean(Up(I).k(:,S,:),3);
        Pknorm = mean(Up(I).Pk(:,S,:)./Up(I).var(:,S,:),3);
        loglog(k.*z_vec(S),Pknorm./z_vec(S),'color',strat_cols(S,:,I),'linewidth',2); hold on;
        set(gca,'fontsize',gca_fontsize)
        if S == 5
            xx = linspace(10^-3,100,100);
            loglog(xx,coeff_AE(1,I)*(xx).^(-1),'k-','linewidth',1.5);
            xx = linspace(0.01,100,100);
            loglog(xx,coeff_ISR(1,I)*(xx).^(-5/3),'k--','linewidth',1.5);
            xline(1,'-','color',[0.7 0.7 0.7],'HandleVisibility','off','linewidth',1.5);
        end
        title({stability_label{I};[]},'FontWeight','normal','FontSize',22);
        ylim([10^-5 200]); xlim([10^-4 200]);
        xticks([10^-4,10^-2,10^0,100]);
    end
    legend(height_labels{S_plot},'location','southwest','interpreter','latex','box','off','fontsize',legend_fontsize);
    if I == 1
        ylabel('$\frac{F_{uu}}{\sigma_{u}^2}$','interpreter','latex','rotation',0,'FontSize',6+label_fontsize);
    end
    text(0.01,0.99,panelLabels{1,I},'units','normalized','interpreter','latex','color','k','fontsize',legend_fontsize,'HorizontalAlignment','left','VerticalAlignment','top');
   

    % --- ROW 2: TEMPERATURE SPECTRA
    subplot(3,3,3+I);
    S_plot = find(Tp(I).Pk(1,:,1));
    for S = S_plot
        k = mean(Tp(I).k(:,S,:),3);
        Pknorm = mean(Tp(I).Pk(:,S,:)./Tp(I).var(:,S,:),3);
        loglog(k.*z_vec(S),Pknorm./z_vec(S),'color',strat_cols(S,:,I),'linewidth',2); hold on;
        if S == 5
            xx = linspace(10^-3,100,100);
            loglog(xx,coeff_AE(1,I)*(xx).^(-1),'k-','linewidth',1.5);
            xx = linspace(0.01,100,100);
            loglog(xx,coeff_ISR(1,I)*(xx).^(-5/3),'k--','linewidth',1.5);
            xline(1,'-','color',[0.7 0.7 0.7],'HandleVisibility','off','linewidth',1.5);
        end
        set(gca,'fontsize',gca_fontsize)
        ylim([10^-5 200]); xlim([10^-4 200]);
        xticks([10^-4,10^-2,10^0,100]);
    end
    legend(height_labels{S_plot},'location','southwest','interpreter','latex','box','off','fontsize',legend_fontsize);
    if I == 1
        ylabel('$\frac{F_{\theta\theta}}{\sigma_{\theta}^2}$','interpreter','latex','rotation',0,'FontSize',6+label_fontsize)
    end
    text(0.01,0.99,panelLabels(2,I),'units','normalized','interpreter','latex','color','k','fontsize',legend_fontsize,'HorizontalAlignment','left','VerticalAlignment','top');
    
    % --- cospectra (ensemble-averaged for plotting)
    subplot(3,3,6+I)
    S_plot = find(UTp(I).Pk(1,:,1));
    us_mean = mean(data(I).ustar);
    Ts_mean = mean(data(I).Tstar);

    for S = S_plot
        kraw = mean(UTp(I).k(:,S,:),3);
        Pkraw = abs(mean(UTp(I).Pk(:,S,:),3));
        Pkrawnorm = abs(mean(UTp(I).Pk(:,S,:)./UTp(I).CuT(:,S,:),3));
        [k,Pk] = smooth_spectra(kraw,Pkraw,500);
        [k,Pknorm] = smooth_spectra(kraw,Pkrawnorm,500);
        loglog(k.*z_vec(S),Pknorm./z_vec(S),'color',strat_cols(S,:,I),'linewidth',2,'HandleVisibility','on'); hold on;
        ind_isr = k.*z_vec(S) < 20 & k.*z_vec(S) > 1;
        [slope(S,I)] = getm(k,Pk,ind_isr);
    end

     if S == 5
            xx = linspace(10^-3,100,100);
            loglog(xx,coeff_AE(1,I)*(xx).^(-1),'k-','linewidth',1.5);
           xx = linspace(0.01,100,100);
            loglog(xx, coeff_ISR(3,I)*(xx).^(-7/3), 'k-.', 'linewidth', 1.5, 'HandleVisibility', 'off');
    xline(1,'-','color',[0.7 0.7 0.7],'HandleVisibility','off','linewidth',1.5);
     end

       legend(height_labels{S_plot},'interpreter','latex','location','southwest','box','off','fontsize',legend_fontsize)
  
    set(gca,'fontsize',gca_fontsize)
    xlabel('$k_xz$','interpreter','latex','FontSize',label_fontsize)
    if I == 1
        ylabel('$\left| \frac{F_{u\theta}}{\overline{ u''\theta_v''}}\right|$','interpreter','latex','rotation',0,'FontSize',6+label_fontsize)
    end
    xline(1,'-','color',[0.7 0.7 0.7],'HandleVisibility','off','linewidth',1.5);
    ylim([10^-5 200]);
    if I == 3, ylim([10^-5 1000]); end
    xlim([10^-4 200]);
    xticks([10^-4,10^-2,10^0,100]);
    text(0.01,0.99,panelLabels(3,I),'units','normalized','interpreter','latex','color','k','fontsize',legend_fontsize,'HorizontalAlignment','left','VerticalAlignment','top');
    
end

exportgraphics(fig, [figures_folder '\Fig_hotwire_spectra_cospectra.pdf'], 'ContentType', 'vector','BackgroundColor','none');

