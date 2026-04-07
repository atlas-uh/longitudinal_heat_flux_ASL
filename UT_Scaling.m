%clear all
clc
CI=3/5;
CR=1.8;
A=(1-CI)/CR;

zeta=-10000:0.001:-0.001;

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

fig = figure(1);
set(gcf,'units','normalized','OuterPosition',[0,0,0.45,0.45]);
co = get(gca,'colororder');
clf
t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile(1)
loglog(-zeta, Phi_m,'-','color',co(1,:),'linewidth',1.5); hold on
loglog(-zeta, Phi_m_KAN,'--','color',co(1,:),'linewidth',1.5)
loglog(-zeta, Phi_h,'-','color',co(2,:),'linewidth',1.5)
loglog(-zeta, Phi_h_KAN,'--','color',co(2,:),'linewidth',1.5)
set(gca,'fontsize',14);
ylabel ('$\phi_m(\zeta),\; \phi_h(\zeta)$','fontweight','bold','FontSize',18,'interpreter','latex')
legend ('$\phi_m$ DDA','$\phi_m$ MOST','$\phi_h$ DDA', '$\phi_h$ MOST','Location','southwest','interpreter','latex')
text(0.01,0.93,'(a)','units','normalized','interpreter','latex','color','k','fontsize',14,'HorizontalAlignment','left','VerticalAlignment','top');
 ylim([10^-1 1])
% xlabel ('$ -\zeta$','fontweight','bold','FontSize',18,'interpreter','latex')
xlim([10^-3 1000])

nexttile(3)
loglog(-zeta, Phi_TKE,'-.','color',co(3,:),'handlevisibility','on','linewidth',1.5)
hold on
loglog(-zeta, Phi_dis,'-','color',co(4,:),'linewidth',1.5)
hold on
loglog(-zeta, Phi_dis_KAN,'--','color',co(4,:),'linewidth',1.5)
hold on
set(gca,'fontsize',14);



ylabel ('$\phi_\varepsilon(\zeta)$, $ \phi_{TKE}(\zeta)$','fontweight','bold','FontSize',18,'interpreter','latex')
legend ('$\phi_{TKE}$','$\phi_{\varepsilon}$ DDA','$\phi_{\varepsilon}$ MOST','Location','best','interpreter','latex')

xlabel ('$ -\zeta$','fontweight','bold','FontSize',18,'interpreter','latex')
text(0.01,0.99,'(b)','units','normalized','interpreter','latex','color','k','fontsize',14,'HorizontalAlignment','left','VerticalAlignment','top');
xlim([10^-3 1000])


hold on

nexttile([2,1])
zz = linspace(0.003,100,10);
loglog(zz, 1.2*zz.^(-2/3),':','color',[0.5 0.5 0.5],'linewidth',2,'HandleVisibility','off')
hold on
loglog(-zeta, R,'k-','linewidth',1.5)
hold on
loglog(-zeta, R_KAN,'k--','linewidth',1.5)
hold on

set(gca,'fontsize',14);
ylabel ('$ R_h=-\overline{u''\theta''}/\overline{w''\theta''}$','fontweight','bold','FontSize',18,'Interpreter','latex')
xlabel ('$-\zeta$','fontweight','bold','FontSize',18,'Interpreter','latex')
text(0.8,1.7,'$1.2 \mathbf{\zeta^{-2/3}}$','fontsize',18,'color',[0.5 0.5 0.5],'Interpreter','latex');
legend ('$R_h$ DDA','$R_h$ MOST','Location','southwest','interpreter','latex')
text(0.01,0.98,'(c)','units','normalized','interpreter','latex','color','k','fontsize',14,'HorizontalAlignment','left','VerticalAlignment','top');

ylim([10^-1 3.2])
xlim([10^-3 1000])
exportgraphics(fig, [figures_folder '\Fig_UT_Scaling.pdf'], 'ContentType', 'vector','BackgroundColor','none');
