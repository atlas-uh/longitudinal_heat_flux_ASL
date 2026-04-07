% LONGITUDINAL HEAT FLUX MAIN PROGRAM 
mkdir('figures');
figures_folder = fullfile(pwd,'figures');

%% --- uT scaling from MOST and DDA 
UT_Scaling; 

%% --- load and analyze sonic data
% --- grass site
load('data_sonic_grass_site.mat');
temp.U = data_sonic_grass_site.u;            
temp.V = data_sonic_grass_site.v;
temp.W = data_sonic_grass_site.w;
[temp] = rotate_sonic(temp);
data_sonic_grass_site.u = temp.U;
data_sonic_grass_site.v = temp.V;
data_sonic_grass_site.w = temp.W;
clear temp
[sonic_grass] = calculate_sonic_stats(data_sonic_grass_site);

% --- SLTEST site 
load('data_sonic_SLTEST.mat');
clear temp
% fill NaNs 
MAX_GAP_SAMPLES = 10;                        
data_sonic_SLTEST.u = fillmissing(data_sonic_SLTEST.u,'linear',1,'EndValues','nearest','MaxGap',MAX_GAP_SAMPLES);
data_sonic_SLTEST.v = fillmissing(data_sonic_SLTEST.v,'linear',1,'EndValues','nearest','MaxGap',MAX_GAP_SAMPLES);
data_sonic_SLTEST.w = fillmissing(data_sonic_SLTEST.w,'linear',1,'EndValues','nearest','MaxGap',MAX_GAP_SAMPLES);
data_sonic_SLTEST.T = fillmissing(data_sonic_SLTEST.T,'linear',1,'EndValues','nearest','MaxGap',MAX_GAP_SAMPLES);

temp.U = data_sonic_SLTEST.u;                
temp.V = data_sonic_SLTEST.v;
temp.W = data_sonic_SLTEST.w;
[temp] = rotate_sonic(temp);

data_sonic_SLTEST.u = temp.U;
data_sonic_SLTEST.v = temp.V;
data_sonic_SLTEST.w = temp.W;

clear temp

[sonic_SLTEST] = calculate_sonic_stats(data_sonic_SLTEST);   

%% --- load and analyze hot-wire data
load('Utah_data.mat') 

for  I = 1:3
    if     I == 1, ii = [605 606];
    elseif I == 2, ii = [581 582 585 586];
    elseif I == 3, ii = [597 598 599];
    end
    data(I).L     = sonic_SLTEST.L(ii);
    data(I).Tstar = sonic_SLTEST.Ts(ii);
    data(I).ustar = sonic_SLTEST.us(ii);
    data(I).wstar = sonic_SLTEST.ws(ii);
end

strat_cols(:,:,2) = flipud(repmat([27,158,119;20 126 96; 15,95,74;11 66 52;8 39 31]./255,2,1));
strat_cols(:,:,1) = flipud(repmat([217,95,2;169 84 2; 126 70 2; 88 54 2; 55 37 2]./255,2,1));
strat_cols(:,:,3) = flipud(repmat([117,112,179;92 89 144;68 67 110;46 46 78;26 26 48]./255,2,1));
fs = 100;
nu = 1.5e-5;
z_vec = [0.0625 0.125 0.25 0.5 1];
height_labels{1} = '0.0625 m'; 
height_labels{2} = '0.125 m'; 
height_labels{3} = '0.25 m'; 
height_labels{4} = '0.5 m'; 
height_labels{5} = '1 m'; 
stability_label{1} = 'unstable';
stability_label{2} = 'near-neutral';
stability_label{3} = 'stable';

% calculate spectra and structure functions
clear Up Tp
nmin = 10;
pow = 13; 
for  I = 1:3
    disp(['--- ' stability_label{I} ' ---'])
    [Up(I),Tp(I)] = process_data_2(data(I).U,data(I).T,pow,fs,data(I).ustar,z_vec,nmin);  
end

% calculate cospectra
clear UTp 
for I = 1:3
     [UTp(I)] = process_cospectra(Up(I),Tp(I),pow,fs); 
end
%% calculate gamma.dUdz and gamma.dTdz from hot-wire data
clear gamma
figure; 
for I = 1:3 % steps through each stability
    % --- gamma_U = dU/dz
    clear dUdz dUdz_2m
    Us = data(I).U;
    
    subplot(2,3,I)
    for k = 1:size(Us,3)
        U = mean(Us(:,:,k));
        z = z_vec;
        p = polyfit(log(z),U,1);
        a = p(1);
        dUdz(:,k) = a./z_vec;   
        dUdz_2m(k) = a./2; % sonic height
        plot(U,z,'ko'); hold on
        zk = linspace(z_vec(1),z_vec(end),20); 
        Uk = polyval(p,log(zk)); 
        plot(Uk,zk,'k-'); xlabel('U'); ylabel('z')
    end 
    gamma(I).dUdz = dUdz;
    gamma(I).dUdz_2m = dUdz_2m;
    % --- gamma T =dT/dz
    Ts = data(I).T;
    clear dTdz dTdz_2m
    subplot(2,3,I+3)
    for k = 1:size(Ts,3)
        T = mean(Ts(:,:,k)); 
        z = z_vec(T>0);
        T = T(T>0);
        p = polyfit(log(z),T,1);
        a = p(1);
        dTdz(:,k) = a./z_vec;  
        dTdz_2m(k) = a./2;
        plot(T,z,'ko'); hold on
        zk = linspace(z_vec(1),z_vec(end),20); 
        Tk = polyval(p,log(zk)); 
        plot(Tk,zk,'k-'); xlabel('T'); ylabel('z')
    end 
    gamma(I).dTdz = dTdz;
    gamma(I).dTdz_2m = dTdz_2m;
end 


%% --- hot-wire: u'T' profile and RuT 
uT_profile;

%% sonic plots 
% --- calculate model Rh
Rh_model;

% --- plot R_uT with stability
sonic_RuT;

% --- plot Rh = uT/wT 
sonic_Rh;

% --- plot spectra
sonic_spectra_cospectra;

% --- C1 (-1 range)
FuT_large_scale_DDA_MOST_coeff; 

%% hot-wire plots

% --- hot-wire: spectra and cospectra 
hotwire_spectra_cospectra_rough; % - case by case to calculate slopes 
hotwire_spectra; % - averaged 

% --- hot-wire: constant scalewise correlation hypothesis 
hotwire_constant_scalewise;

% --- hot-wire: cospectral budget
cospectral_budget_solution;

%% --- appendix: hotwire and sonic comparison
hotwire_sonic_comparison;
