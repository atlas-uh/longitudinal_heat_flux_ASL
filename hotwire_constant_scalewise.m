fig = figure(9); clf; 
set(gcf,'units','normalized','OuterPosition',[0,0,0.6,0.5])
beta_cutoff = [2.5 1 0.125]; 
deltas = [350 60 30]; % boundary layer heights
tl = tiledlayout(2,3);
tl.TileSpacing = "compact";
tl.Padding     = "compact";
mult = 10; 
for I = 1:3
    Cs = UTp(I).Cs; 
    S_plot = find(UTp(I).Pf(1,:,1));

    % -----------
    nexttile(tl,I);
    for S = S_plot
        [k_smooth,C_smooth]=smooth_spectra(mean(UTp(I).k(:,S,:),3),mean(Cs(:,S,:),3),500);
        k = k_smooth;
        C = C_smooth;
        loglog(k.*z_vec(S),C,'linewidth',1.5,'color',strat_cols(S,:,I)); hold on;
        k1 = 1/(mult*mean(Up(I).integral_length_scale_auto(:,S,:),3)); 
        k2 = 1/((1/mult)*mean(Up(I).integral_length_scale_auto(:,S,:),3)); 
        k_ISR_end = f_cutoff.*(2*pi)./mean(Up(I).mean(:,:,:),3);
        k3 = k_ISR_end(S);
        ind = find(k>k2 & k<k3);
        [slope] = get_spectral_exponent(k,C,ind); 
        m(I).C.ISR(S) = slope;
        ind = find(k<k2);
        [slope] = get_spectral_exponent(k,C,ind); 
        m(I).C.AE(S) = slope;
        ylim([0.07 1]);
   end
    xline(1,'-','color',[0.7 0.7 0.7],'HandleVisibility','off','linewidth',2);
    zL = z_vec./mean(data(I).L);
     %xlim([10^-4 200])
    
   if I == 2 
        set(gca,'fontsize',gca_fontsize,'XScale', 'log', 'XTick', logspace(-4, 2, 4)); 
   else
           set(gca,'fontsize',gca_fontsize,'XScale', 'log', 'XTick', logspace(-4, 2, 7)); 

   end
   %     leg = legend(num2str(zL(1),'%.3f'),num2str(zL(2),'%.3f'),num2str(zL(3),'%.3f'),num2str(zL(4),'%.3f'),num2str(zL(5),'%.3f'),'location','southwest','fontsize',12);
   %     title(leg,'z/L')


    leg = legend(num2str(z_vec(S_plot)','%.3f'),'location','southwest','fontsize',12);
    title(leg,'$z$ (m)','interpreter','latex');

    
    xlabel('$k_xz$','interpreter','latex','fontsize',20);
 title({stability_label{I};[]},'FontWeight','normal','FontSize',18);


    nexttile(tl,I+3)
   for S = S_plot
        [k_smooth,C_smooth]=smooth_spectra(mean(UTp(I).k(:,S,:),3),mean(Cs(:,S,:),3),500);
        k = k_smooth;
        C = C_smooth;
        d = 0.0424;% sensor separation distance
        loglog(k.*d,C,'linewidth',1.5,'color',strat_cols(S,:,I)); hold on;
        k1 = 1/(mult*mean(Up(I).integral_length_scale_auto(:,S,:),3)); 
        k2 = 1/((1/mult)*mean(Up(I).integral_length_scale_auto(:,S,:),3)); 
        k_ISR_end = f_cutoff.*(2*pi)./mean(Up(I).mean(:,:,:),3);
        k3 = k_ISR_end(S);
        ind = find(k>k2 & k<k3);
        [slope] = get_spectral_exponent(k,C,ind); 
        m(I).C.ISR(S) = slope;
        ind = find(k<k2);
      
        [slope] = get_spectral_exponent(k,C,ind); 
        m(I).C.AE(S) = slope;
        ylim([0.07 1]);

    end
    xline(1,'-','color',[0.7 0.7 0.7],'HandleVisibility','off','linewidth',2);
    zL = z_vec./mean(data(I).L);
    set(gca,'fontsize',gca_fontsize,'XScale', 'log', 'XTick', logspace(-4, 2, 7)); 
    xlabel('$k_xd$','interpreter','latex','fontsize',20);
end
 
  ylabel(tl,'${F_{u \theta}}/{\left(F_{uu}^{1/2}F_{\theta \theta}^{1/2}\right)}$','interpreter','latex','fontsize',label_fontsize);



exportgraphics(fig, [figures_folder '\Fig_CSC.pdf'], 'ContentType', 'vector','BackgroundColor','none');
