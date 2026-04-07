
%% --- store all data points across runs for each stability case 
clear temp
for I = 1:3 
    % clear variables
    uT_all=[]; z_all=[]; uT_norm_all = [];
    clear uT uT_norm 
    for E = 1:size(Up(I).fluc,3) % for each run 
        % --- hot-wire
        uprime = Up(I).fluc(:,:,E);
        Tprime = Tp(I).fluc(:,:,E);
        uT = mean(uprime.*Tprime); 
        uT_norm = mean(uprime.*Tprime)./sqrt(Up(I).var(:,:,E))./sqrt(Tp(I).var(:,:,E));     

        % --- sonic
        if     I == 1, ii = [605 606];
        elseif I == 2, ii = [581 582 585 586];
        elseif I == 3, ii = [597 598 599];
        end
        uT_2m = [sonic_SLTEST.CuT(ii(E))];
        sigu_2m = [sqrt(sonic_SLTEST.Cuu(ii(E)))];
        sigT_2m = [sqrt(sonic_SLTEST.CTT(ii(E)))];
        uT(6) = uT_2m; 
        uT_norm(6) = uT_2m./sigT_2m./sigu_2m; 
        z = [z_vec 2]; 

        uT_norm_all = [uT_norm_all uT_norm(:)];
        uT_all = [uT_all uT(:)];
        z_all = [z_all z(:)];

    end
    temp(I).uT_norm = uT_norm_all; 
    temp(I).uT = uT_all; 
    temp(I).z = z_all; 
end


%% --- plot
fig = figure (3); clf; 
set(gcf,'units','normalized','OuterPosition',[0,0,0.55,0.4]);

CI = 3/5; CR = 1.8; 
alpha = 0.5;
legend_fontsize = 14;
gca_fontsize = 16;
label_fontsize = 20; 
mstyle = ['v';'o';'s'];
   
% --- (a)
subplot(1,2,1)
for I = 1:3 
    uT_norm_all = temp(I).uT_norm; uT_norm = uT_norm_all(~isnan(uT_norm_all(:,1)),:);
    z = temp(I).z; z = z(~isnan(uT_norm_all(:,1)),:);
    
    % --- plot each case as thin lines 
    for E = 1:size(z,2)
          semilogy(uT_norm(:,E),z(:,E),'-','color',[strat_cols(4,:,I) 0.3],'HandleVisibility','off','linewidth',1.5); hold on
  
    end
    x = median(uT_norm,2);
    y = mean(z,2);
    semilogy(x,y,mstyle(I,:),'color',strat_cols(4,:,I),'markerfacecolor',strat_cols(4,:,I),'markersize',8,'HandleVisibility','on','linewidth',1.5); hold on
    % calculate p value (if >0.05, then slope is not statistically significantly different than 0) 
    mdl = fitlm(y, x);
    pVals_norm(I) = mdl.Coefficients.pValue(2); 


end 
% --- (b)
subplot(1,2,2); 

for I = 1:3
    uT_all = temp(I).uT; uT = uT_all(~isnan(uT_all(:,1))& abs(uT_all(:,1))>0,:);
    z = temp(I).z; z = z(~isnan(uT_all(:,1))& abs(uT_all(:,1))>0,:);
    
 
    % plot each case lightly
    for E = 1:size(z,2)
          semilogy(uT(:,E),z(:,E),'-','color',[strat_cols(4,:,I) 0.3],'HandleVisibility','off','linewidth',1.5); hold on;
    end
 
    % plot median (representative)
    x = median(uT,2);
    y = mean(z,2);
    semilogy(x,y,mstyle(I,:),'color',strat_cols(4,:,I),'markerfacecolor',strat_cols(4,:,I),'markersize',8,'HandleVisibility','on','linewidth',1.5); hold on
    % calculate p value (if >0.05, then slope is not statistically
    % significantly different than 0) 
    mdl = fitlm(y, x);
    pVals(I) = mdl.Coefficients.pValue(2); 
end

% --- plot formatting 
subplot(1,2,1)
xline([0 0],'k-','handlevisibility','off'); hold on;
text(0.01,0.99,'(a)','units','normalized','interpreter','latex','color','k','fontsize',legend_fontsize+4,'HorizontalAlignment','left','VerticalAlignment','top');
set(gca,'fontsize',gca_fontsize);
xlabel('$R_{u\theta}=\overline{ u''\theta_v''}/\left( \sigma_u \sigma_\theta\right)$','interpreter','latex','fontsize',label_fontsize);
ylabel({'$z\;\;$   '},'interpreter','latex','fontsize',label_fontsize+4,'rotation',0);%plot(uT_2m,2,'o'); 
xlim([-1 1]);
ylim([0.05 3]);
leg_cell = [];
leg_cell{1} = stability_label{1};
leg_cell{2} = stability_label{2};
leg_cell{3} = stability_label{3};
lgd = legend(leg_cell,'fontsize',legend_fontsize,'numcolumns',1,'orientation','horizontal','location','northeast','box','on');
pos = lgd.Position; 
lgd.Position = [0.2209 0.3590 pos(3) pos(4)];


subplot(1,2,2)
xlim([-0.3 0.3]);
ylim([0.05 3]);
set(gca,'fontsize',gca_fontsize);
xlabel('$\overline{ u''\theta_v'' }$','interpreter','latex','fontsize',label_fontsize);
xline([0 0],'k-','handlevisibility','off'); hold on;
leg_cell = [];
text(0.01,0.99,'(b)','units','normalized','interpreter','latex','color','k','fontsize',legend_fontsize+4,'HorizontalAlignment','left','VerticalAlignment','top');


exportgraphics(fig, [figures_folder '\Fig_uT_prof.pdf'], 'ContentType', 'vector','BackgroundColor','none');


%% --- report z-variation statistics
stability_label_short = {'Unstable', 'Near-neutral', 'Stable'};

fprintf('\n--- Linear regression of median profile vs log(z) ---\n');
fprintf('p > 0.05 means slope is NOT significantly different from zero\n\n');

fprintf('%-15s %-25s %-25s\n', 'Stability', 'R_utheta (norm)', 'u''theta'' (raw)');
fprintf('%s\n', repmat('-',1,65));
for I = 1:3
    fprintf('%-15s p = %-20.3f p = %-20.3f\n', ...
        stability_label_short{I}, pVals_norm(I), pVals(I));
end
