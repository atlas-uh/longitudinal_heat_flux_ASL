% --- Define base colors
blue  = [0.00 0.30 0.55];
orange = [0.70 0.25 0.05];
label_fontsize = 20;
gca_fontsize = 15;
legend_fontsize = 15;
m_size = 8; % marker size

clear eps;
Rh_SLTEST = -sonic_SLTEST.CuT./(10*eps+sonic_SLTEST.CwT);
Rh_grass = -sonic_grass.CuT./(10*eps+sonic_grass.CwT);

fig = figure(5); close(fig); 
fig = figure(5); 
co = get(gca,'ColorOrder');
set(gcf,'units','normalized','OuterPosition',[0,0,0.64,1])
tl = tiledlayout(3,2);
tl.TileSpacing = "compact";
alpha = 0.4;
scatter_size = 20; 

%% (a) Rh vs zeta 

subplot(3,2,1)

xlimits = [-1 7];
ylimits = [-5 20];

hold on
xline(0,'k-','HandleVisibility','off');
yline(0,'k-','HandleVisibility','off');


alpha1 = 0.2;              % SLTEST transparency (smaller = lighter)
alpha2 = 0.30;              % Grass transparency
scatter_size = 20;
SLTEST_marker = 'o';
grass_marker = 's';

x1 = -sonic_SLTEST.zeta(:);
y1 = Rh_SLTEST(:);
tol = 1e-10;

x2 = -sonic_grass.zeta(:);
y2 = Rh_grass(:);
mask1 = ~(abs(x1) < tol & abs(y1) < tol);
mask2 = ~(abs(x2) < tol & abs(y2) < tol);

% --- SLTEST (blue)
s1 = scatter(x1(mask1), y1(mask1), scatter_size, 'filled', 'MarkerFaceColor', blue, 'MarkerEdgeColor','none');
s1.MarkerFaceAlpha = alpha1;  
s1.Marker = SLTEST_marker; 
% --- Grass (orange)
s2 = scatter(x2(mask2), y2(mask2), scatter_size, 'filled','MarkerFaceColor', orange, 'MarkerEdgeColor', 'none');
s2.MarkerFaceAlpha = alpha2;
s2.Marker = grass_marker;




% --- Formatting
xlim(xlimits);
ylim(ylimits);
set(gca,'FontSize',gca_fontsize);

xlabel('$-\zeta$','interpreter','latex','FontSize',label_fontsize);
ylabel('$R_h=-\overline{u''\theta_v''}/\overline{w''\theta_v''}$', ...
       'interpreter','latex','Rotation',90,'FontSize',label_fontsize);

box on
text(0.01,0.99,'(a)','units','normalized','interpreter','latex', ...
     'color','k','fontsize',legend_fontsize+2, ...
     'HorizontalAlignment','left','VerticalAlignment','top');
lgd = legend('SLTEST','Grass Clearing','location','southeast', ...
    'box','on','fontsize',legend_fontsize,'numColumns',2,'interpreter','latex');

% --- inset
axes('Position',[.22 .81 .23 .1])
box on

% --- SLTEST (blue)
s1 = scatter(x1(mask1), y1(mask1), scatter_size, 'filled', 'MarkerFaceColor', blue, 'MarkerEdgeColor','none');
s1.MarkerFaceAlpha = alpha1;   
s1.Marker = SLTEST_marker;   
hold on

% --- Grass (orange)
s2 = scatter(x2(mask2), y2(mask2), scatter_size, 'filled','MarkerFaceColor', orange, 'MarkerEdgeColor', 'none');
s2.MarkerFaceAlpha = alpha2;
s2.Marker = grass_marker;


xlim([-0.5 0.5]);
ylim([-1 6])
xline([0 0],'k-');
yline([0 0],'k-');
box on
set(gca,'fontsize',gca_fontsize)




%% (b) uT vs wT for near-neutral

subplot(3,2,2)

% select NN only 
ii_SLTEST = find(abs(sonic_SLTEST.zeta)<0.05); 
ii_grass = find(abs(sonic_grass.zeta)<0.05); 


xlimits = [-0.05 0.3];
ylimits = [-1 0.5];

hold on
xline(0,'k-','HandleVisibility','off');
yline(0,'k-','HandleVisibility','off');

x1 = sonic_SLTEST.CwT(ii_SLTEST)';
y1 = sonic_SLTEST.CuT(ii_SLTEST)';

x2 = sonic_grass.CwT(ii_grass)';
y2 = sonic_grass.CuT(ii_grass)';

mask1 = ~(abs(x1) < tol & abs(y1) < tol);
mask2 = ~(abs(x2) < tol & abs(y2) < tol);


%plot fit line first so it's behind data
xfit = linspace(xlimits(1),xlimits(2),100);
p_1 = polyfit(x1,y1,1);
plot(xfit,polyval(p_1,xfit),'-','color',blue,'linewidth',1); hold on; 

p_2 = polyfit(x2,y2,1);
plot(xfit,polyval(p_2,xfit),'-','color',orange,'linewidth',1);

% -3 refernce 
plot(xfit,-3*(xfit),'--','color','k','linewidth',1);


% --- SLTEST (blue)
s1 = scatter(x1(mask1), y1(mask1), scatter_size, 'filled', 'MarkerFaceColor', blue, 'MarkerEdgeColor','none');
s1.MarkerFaceAlpha = alpha1;    
s1.Marker = SLTEST_marker;   

% --- Grass (orange)
s2 = scatter(x2(mask2), y2(mask2), scatter_size, 'filled','MarkerFaceColor', orange, 'MarkerEdgeColor', 'none');
s2.MarkerFaceAlpha = alpha2;
s2.Marker = grass_marker;


% --- Formatting
xlim(xlimits);
ylim(ylimits);
set(gca,'FontSize',gca_fontsize);

xlabel('$-\zeta$','interpreter','latex','FontSize',label_fontsize);
ylabel('$R_h=-\overline{u''\theta_v''}/\overline{w''\theta_v''}$', ...
       'interpreter','latex','Rotation',90,'FontSize',label_fontsize);

ylimits = [-1 0.5];
plot(xlimits,[0 0],'k-'); hold on;
plot([0 0],ylimits,'k-'); 
xlim(xlimits);
ylim(ylimits);
set(gca,'FontSize',gca_fontsize);


xlabel('$\overline{ w''\theta_v'' }$','interpreter','latex','FontSize',label_fontsize,'rotation',0);
ylabel('$\overline{ u''\theta_v'' }$','interpreter','latex','fontsize',label_fontsize);
text(0.01,0.99,'(b) for $|\zeta|<0.05$','units','normalized','interpreter','latex','color','k','fontsize',legend_fontsize+2,'HorizontalAlignment','left','VerticalAlignment','top');

legend(['SLTEST $R_h=' num2str(-p_1(1),'%.2f') '$'],...
       ['Grass $R_h=' num2str(-p_2(1),'%.2f') '$'],...
       'interpreter','latex','location','southwest','fontsize',gca_fontsize-1,'box','on');
box on
text(0.26,-0.655,'$-3$','interpreter','latex', ...
     'color','k','fontsize',legend_fontsize+2, ...
     'HorizontalAlignment','left','VerticalAlignment','top');




% --- inset 
axes('Position',[.8,.82,.09 .09])
%plot fit line first so it's behind data
xfit = linspace(xlimits(1),xlimits(2),100);
p_1 = polyfit(x1,y1,1);
plot(xfit,polyval(p_1,xfit),'-','color',blue,'linewidth',1); hold on; 

p_2 = polyfit(x2,y2,1);
plot(xfit,polyval(p_2,xfit),'-','color',orange,'linewidth',1);

plot(xfit,-3*(xfit),'--','color','k','linewidth',1);

% --- SLTEST (blue)
s1 = scatter(x1(mask1), y1(mask1), scatter_size, 'filled', 'MarkerFaceColor', blue, 'MarkerEdgeColor','none');
s1.MarkerFaceAlpha = alpha1;   
s1.Marker = SLTEST_marker;  
hold on;

% --- Grass (orange)
s2 = scatter(x2(mask2), y2(mask2), scatter_size, 'filled','MarkerFaceColor', orange, 'MarkerEdgeColor', 'none');
s2.MarkerFaceAlpha = alpha2;
s2.Marker = grass_marker;

xlim([-0.02 0.02]);
ylim([-0.1 0.1])
yticks([-0.1 0 0.1])
xline([0 0],'k-');
yline([0 0],'k-');
set(gca,'fontsize',gca_fontsize)


%% (c) Rh /(u*^2/w*^2)
subplot(3,2,3)

xlimits = [-1 7];
ylimits = [-5 2];

hold on
xline(0,'k-','HandleVisibility','off');
yline(0,'k-','HandleVisibility','off');

us_ws_SLTEST = (sonic_SLTEST.us.^2)./(sonic_SLTEST.ws.^2); 
Rh_DDA_SLTEST = Rh_SLTEST./us_ws_SLTEST;
us_ws_grass = (sonic_grass.us.^2)./(sonic_grass.ws.^2);
Rh_DDA_grass = Rh_grass./us_ws_grass;


x1 = -sonic_SLTEST.zeta(:);
y1 = real(Rh_DDA_SLTEST)';

x2 = -sonic_grass.zeta(:);
y2 = real(Rh_DDA_grass)';

mask1 = ~(abs(x1) < tol & abs(y1) < tol);
mask2 = ~(abs(x2) < tol & abs(y2) < tol);


% --- SLTEST (blue)
s1 = scatter(x1(mask1), y1(mask1), scatter_size, 'filled', 'MarkerFaceColor', blue, 'MarkerEdgeColor','none');
s1.MarkerFaceAlpha = alpha1;    
s1.Marker = SLTEST_marker;    
hold on;

% --- Grass (orange)
s2 = scatter(x2(mask2), y2(mask2), scatter_size, 'filled','MarkerFaceColor', orange, 'MarkerEdgeColor', 'none');
s2.MarkerFaceAlpha = alpha2;
s2.Marker = grass_marker;

% --- Formatting
xlim(xlimits);
ylim(ylimits);
set(gca,'FontSize',gca_fontsize);

ylabel('$R_h/\left(u_*^2/w_*^2\right)$','interpreter','latex','fontsize',label_fontsize,'rotation',90);
xlabel('$-\zeta$','interpreter','latex','FontSize',label_fontsize);

box on
text(0.01,0.99,'(c)','units','normalized','interpreter','latex', ...
     'color','k','fontsize',legend_fontsize+2, ...
     'HorizontalAlignment','left','VerticalAlignment','top');


% --- inset
axes('Position',[.205 .445 .25 .09])

% --- SLTEST (blue)
s1 = scatter(x1(mask1), y1(mask1), scatter_size, 'filled', 'MarkerFaceColor', blue, 'MarkerEdgeColor','none');
s1.MarkerFaceAlpha = alpha1;    
s1.Marker = SLTEST_marker;   
hold on;

% --- Grass (orange)
s2 = scatter(x2(mask2), y2(mask2), scatter_size, 'filled','MarkerFaceColor', orange, 'MarkerEdgeColor', 'none');
s2.MarkerFaceAlpha = alpha2;
s2.Marker = grass_marker;

xlim([-0.5 2]);
ylim([-1 2])
xline([0 0],'k-');
yline([0 0],'k-');
set(gca,'fontsize',gca_fontsize)
box on
yticks([-2  0  2])



% plateau 

% Pool both datasets (use inset range only)
x_all = [x1(:); x2(:)];
y_all = [y1(:); y2(:)];

% x-range
xmin = 0.25;    
xmax = 20;
mr = (x_all >= xmin) & (x_all <= xmax);
xf = x_all(mr);
yf = y_all(mr);

plateau = median(yf) ;

hold on
yline(plateau, 'k:', 'LineWidth', 2, ...
      'DisplayName', sprintf('asymptote \\approx %.2f', plateau));

text(-0.5,plateau,[num2str(plateau,'%.2f') ' '],'interpreter','latex','color','k',...
    'fontsize',legend_fontsize,'HorizontalAlignment','right','VerticalAlignment','middle');



%% calculate collapse 
xi0 = 0.05;

% near-neutral masks (use zeta, not -zeta; abs() makes sign irrelevant)
m1 = abs(sonic_SLTEST.zeta(:)) < xi0;
m2 = abs(sonic_grass.zeta(:)) < xi0;

% --- (a) MOST scaling ordinate: Rh
yA1 = Rh_SLTEST(:);
yA2 = Rh_grass(:);

% --- (c) DDA scaling ordinate: Rh/(u*^2/w*^2)
yC1 = real(Rh_DDA_SLTEST(:));
yC2 = real(Rh_DDA_grass(:));

% pooled vectors (both datasets together)
yA = [yA1(m1); yA2(m2)];
yC = [yC1(m1); yC2(m2)];

yA = yA(~isnan(yA));
yC = yC(~isnan(yC));

% --- Quartiles
Q1_A = prctile(yA,25);   Q3_A = prctile(yA,75);
Q1_C = prctile(yC,25);   Q3_C = prctile(yC,75);

% --- IQR (absolute spread of the central 50%)
IQR_A = Q3_A - Q1_A;
IQR_C = Q3_C - Q1_C;

% --- Normalized IQR (scale-fair; easiest to explain)
IQRn_A = IQR_A / (abs(median(yA)) + eps);
IQRn_C = IQR_C / (abs(median(yC)) + eps);

% --- CQV (alternative normalized robust spread)
CQV_A = IQR_A / (abs(Q3_A + Q1_A) + eps);
CQV_C = IQR_C / (abs(Q3_C + Q1_C) + eps);

% --- Improvements (% reduction from (a) -> (c))
improve_IQR  = 100*(1 - IQR_C/IQR_A);
improve_IQRn = 100*(1 - IQRn_C/IQRn_A);
improve_CQV  = 100*(1 - CQV_C/CQV_A);


%% (d) Rh with predicitons
%nexttile(tl,4)
subplot(3,2,4)

xlimits = [10^-3 20];
ylimits = [0.001 50];
hold on

x1 = -sonic_SLTEST.zeta(:);
y1 = real(Rh_SLTEST)';

x2 = -sonic_grass.zeta(:);
y2 = real(Rh_grass)';

mask1 = ~(abs(x1) < tol & abs(y1) < tol);
mask2 = ~(abs(x2) < tol & abs(y2) < tol);

% --- plot predicted 
CI=3/5;
CR=1.8;
A=(1-CI)/CR;

zeta=-100:0.00001:-0.0001;

sig_u=2.7*(1-0*3*zeta).^(1/3);
sig_v=2.4*(1-0*3*zeta).^(1/3);
sig_w=1.3*(1-1*3*zeta).^(1/3);

Phi_TKE=0.5*(sig_u.^2+sig_v.^2+sig_w.^2);

Phi_m_KAN=(1+16*(-zeta)).^(-1/4);
Phi_h_KAN=(1+16*(-zeta)).^(-1/2);
Phi_dis_KAN=(1+0.5*(-zeta).^(3/2)).^(2/3);

Phi_m=((1+0.625*(-zeta).^2)./(1+7.5*(-zeta))).^(1/3);
Phi_h=0.64*((3+2.5*(-zeta))./(1+10*(-zeta)+50*(-zeta).^2)).^(1/3);
Phi_dis=0.4*(10+7.5*(-zeta)+6.25*(-zeta).^2)./(4+2.5*(-zeta));

R_KAN=A*(Phi_TKE.*Phi_m_KAN).*(1+Phi_h_KAN./Phi_m_KAN)./Phi_dis_KAN;
R=A*(Phi_TKE.*Phi_m).*(1+Phi_h./Phi_m)./Phi_dis;


loglog(-zeta,R,'k-','LineWidth',1); hold on;
loglog(-zeta,R_KAN,'k--','LineWidth',1)
legend('DDA','MOST','location','southwest');
box on
% --- Scatter

% --- SLTEST (blue)
s1 = scatter(x1(mask1), y1(mask1), scatter_size, 'filled', 'MarkerFaceColor', blue, 'MarkerEdgeColor','none','HandleVisibility','off');
s1.MarkerFaceAlpha = alpha1;   
s1.Marker = SLTEST_marker;    
hold on;

% --- Grass (orange)
s2 = scatter(x2(mask2), y2(mask2), scatter_size, 'filled','MarkerFaceColor', orange, 'MarkerEdgeColor', 'none','HandleVisibility','off');
s2.MarkerFaceAlpha = alpha2;
s2.Marker = grass_marker;


xl = [0.001 100];
loglog(xl,1.2*xl.^(-2/3),'k:','LineWidth',2,'HandleVisibility','off')
set(gca,'fontsize',gca_fontsize)
ylabel('$R_h$','interpreter','latex','fontsize',label_fontsize,'rotation',90);
xlabel('$-\zeta$','interpreter','latex','FontSize',label_fontsize);
text(0.01,0.99,'(d)','units','normalized','interpreter','latex','color','k','fontsize',legend_fontsize+2,'HorizontalAlignment','left','VerticalAlignment','top');
set(gca, 'XScale', 'log', 'YScale', 'log');
text(0.06,15,'$-2/3$','interpreter','latex','color','k','fontsize',legend_fontsize);
xlim(xlimits)
ylim(ylimits)



%% (e) stable scaling

subplot(3,2,5.5)

xlimits = [-1 7];
ylimits = [-0.5 2];
hold on
xline(0,'k-','HandleVisibility','off');
yline(0,'k-','HandleVisibility','off');


DO_scale_SLTEST = sonic_SLTEST.eps_son.D2(:)./sonic_SLTEST.beta(:);
uT_DO_SLTEST = sonic_SLTEST.CuT(:)./DO_scale_SLTEST(:); 

DO_scale_grass = sonic_grass.eps_son.D2(:)./sonic_grass.beta(:);
uT_DO_grass = sonic_grass.CuT(:)./DO_scale_grass(:); 


x1 = -sonic_SLTEST.zeta(:);
y1 = uT_DO_SLTEST(:);

x2 = -sonic_grass.zeta(:);
y2 = uT_DO_grass(:);

mask1 = ~(abs(x1) < tol & abs(y1) < tol);
mask2 = ~(abs(x2) < tol & abs(y2) < tol);

% --- Scatter

% --- SLTEST (blue)
s1 = scatter(x1(mask1), y1(mask1), scatter_size, 'filled', 'MarkerFaceColor', blue, 'MarkerEdgeColor','none','HandleVisibility','off');
s1.MarkerFaceAlpha = alpha1;     
s1.Marker = SLTEST_marker;    
hold on;

% --- Grass (orange)
s2 = scatter(x2(mask2), y2(mask2), scatter_size, 'filled','MarkerFaceColor', orange, 'MarkerEdgeColor', 'none','HandleVisibility','off');
s2.MarkerFaceAlpha = alpha2;
s2.Marker = grass_marker;


% --- Formatting
xlim([-1.0 0]);
ylim([-0.25 2.4]);
set(gca,'fontsize',gca_fontsize, 'box','on')
xlabel('$-\zeta$','interpreter','latex','FontSize',label_fontsize);
ylabel('$\overline{u''\theta_v''}/\left(U_{DO}\theta_{DO}\right)$','interpreter','latex','Rotation',90,'FontSize',label_fontsize);
text(0.01,0.99,'(e)','units','normalized','interpreter','latex','color','k','fontsize',legend_fontsize+2,'HorizontalAlignment','left','VerticalAlignment','top');


% --- inset 
axes('Position',[.41 .225 .25 .09])
% --- Scatter

% --- SLTEST (blue)
s1 = scatter(x1(mask1), y1(mask1), scatter_size, 'filled', 'MarkerFaceColor', blue, 'MarkerEdgeColor','none','HandleVisibility','off');
s1.MarkerFaceAlpha = alpha1;   
s1.Marker = SLTEST_marker;   
hold on;

% --- Grass (orange)
s2 = scatter(x2(mask2), y2(mask2), scatter_size, 'filled','MarkerFaceColor', orange, 'MarkerEdgeColor', 'none','HandleVisibility','off');
s2.MarkerFaceAlpha = alpha2;
s2.Marker = grass_marker;


xlim([-0.2 0]);
yticks(0:0.2:0.4)
ylim([0.0 0.4]);
xline([0 0],'k-');
yline([0 0],'k-');
yline(-0.5,'k--');
yline(-3,'k--');
box on
set(gca,'fontsize',gca_fontsize)

% --- binned medians with IQR error bars on inset
edges_e = -1:0.02: 0; 
bin_centers_e = (edges_e(1:end-1) + edges_e(2:end)) / 2;

med_e   = NaN(length(bin_centers_e), 1);
q1_e    = NaN(length(bin_centers_e), 1);
q3_e    = NaN(length(bin_centers_e), 1);

x_pool = [-sonic_SLTEST.zeta(:); -sonic_grass.zeta(:)];
y_pool = [uT_DO_SLTEST(:); uT_DO_grass(:)];
valid_pool = isfinite(x_pool) & isfinite(y_pool) & abs(y_pool) > 0;

for b = 1:length(bin_centers_e)
    in_bin = valid_pool & x_pool >= edges_e(b) & x_pool < edges_e(b+1);
    vals = y_pool(in_bin);
    if length(vals) >= 3
        med_e(b) = median(vals);
        q1_e(b)  = prctile(vals, 25);
        q3_e(b)  = prctile(vals, 75);
    end
end

ok_e = ~isnan(med_e);
errorbar(bin_centers_e(ok_e), med_e(ok_e), ...
    med_e(ok_e) - q1_e(ok_e), ...
    q3_e(ok_e) - med_e(ok_e), ...
    'k-o', 'linewidth', 1, 'MarkerSize',0.2 , 'MarkerFaceColor', 'k', 'CapSize', 4)


exportgraphics(fig, [figures_folder '\Fig_Rh.pdf'], 'ContentType', 'vector','BackgroundColor','none');


