%% 
kappa = 0.4; 

% calculate phi_TKE(0) from SLTEST sonic
y = sonic_SLTEST.TKE./(sonic_SLTEST.us.^2);
x = - sonic_SLTEST.zeta; 
mask =  abs(x)< 0.01 & ~isnan(x) & ~isnan(y) & y<50; %
xx = x(mask); 
yy = y(mask);

phi_TKE_0 = mean(yy);

% calcualte phi_m(0), phi_h(0), and phi_eps(0) from near-neutral set (I=2) 

ii = [581,582,585,586];
I = 2;
clear phi_h phi_m phi_eps phi_TKE
for k = 1:length(ii)
    phi_h(k)   = (kappa*1) .* gamma(I).dTdz(5,k) ./ sonic_SLTEST.Ts(ii(k));
    phi_m(k)   = (kappa*1) .* gamma(I).dUdz(5,k) ./ sonic_SLTEST.us(ii(k));
    eps = squeeze(Up(I).eps.D3(1,5,:));
    phi_eps(k) = eps(k) .* kappa .* 1 ./ (sonic_SLTEST.us(ii(k)).^3);
end

Rh_MOST  = -((1-(3/5))/1.8)*(6.7)*(1+1);
Rh_model_SLTEST = mean(-((1-(3/5))/1.8)*((phi_TKE_0).*(phi_m)./(phi_eps)).*(1+((phi_h)./(phi_m))));

fprintf('phi_TKE_0 = %.2f\n', phi_TKE_0)
fprintf('phi_m   = %.2f\n', mean(phi_m))
fprintf('phi_h   = %.2f\n', mean(phi_h))
fprintf('phi_eps = %.2f\n', mean(phi_eps))
fprintf('R_h = %.2f\n', Rh_model_SLTEST)

