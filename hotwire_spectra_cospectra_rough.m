clear m

for I = 1:3
    figure;

    for E = 1:size(data(I).U,3)
        for S = 1:5

            % --- band limits 
            Lint = Up(I).integral_length_scale_auto(1,S,E);
            Ubar = Up(I).mean(1,S,E);
            
            mult = 10;
            k_ISR_start = mult / Lint;
            f_cutoff = 20;
            k_ISR_end   = f_cutoff * (2*pi) / Ubar;


            % -------------------------
            % Row 1: uu spectrum
            % -------------------------
            subplot(3,5,S);

            k    = Up(I).k(:,S,E);
            P    = Up(I).Pk(:,S,E);
            Puu  = P;
          
            loglog(k,P,'k'); hold on;

            ind = (k > k_ISR_start) & (k < k_ISR_end);
            if any(ind)
                loglog(k(ind),P(ind),'r'); hold on;
                m(I).phi_uu.ISR(E,S) = get_spectral_exponent(k,P,ind);
            else
                m(I).phi_uu.ISR(E,S) = NaN;
            end

            ind_AE = (k < k_ISR_start);
            if any(ind_AE)
                m(I).phi_uu.AE(E,S) = get_spectral_exponent(k,P,ind_AE);
            else
                m(I).phi_uu.AE(E,S) = NaN;
            end

            if E == 1
                loglog(k,0.1*k.^(-1),'k-');
                loglog(k,0.1*k.^(-5/3),'b--','linewidth',1.5);
            end

            % -------------------------
            % Row 2: TT spectrum
            % -------------------------
            subplot(3,5,5+S);

            k    = Tp(I).k(:,S,E);
            P    = Tp(I).Pk(:,S,E);
            PTT  = P; 


            loglog(k,P,'k'); hold on;

            ind = (k > k_ISR_start) & (k < k_ISR_end);
            if any(ind)
                loglog(k(ind),P(ind),'r'); hold on;
                m(I).phi_TT.ISR(E,S) = get_spectral_exponent(k,P,ind);
            else
                m(I).phi_TT.ISR(E,S) = NaN;
            end

            ind_AE = (k < k_ISR_start);
            if any(ind_AE)
                m(I).phi_TT.AE(E,S) = get_spectral_exponent(k,P,ind_AE);
            else
                m(I).phi_TT.AE(E,S) = NaN;
            end

            if E == 1
                loglog(k,0.1*k.^(-1),'k-');
                loglog(k,0.1*k.^(-5/3),'b--','linewidth',1.5);
            end

            % -------------------------
            % Row 3: uT cospectrum
            % -------------------------
            subplot(3,5,10+S);

            
            TbarK = Tp(I).mean(1,S,E);
            beta  = 9.81 / TbarK;

            NBV(S)   = sqrt(beta * gamma(I).dTdz(S,E));
            k_NBV(S) = NBV(S) * (2*pi) / Ubar;

            % cospectrum k must be the column for this sensor
            k    = UTp(I).k(:,S,E);
            Praw    = abs(UTp(I).Pk(:,S,E));

            [k_smooth, P_smooth] = smooth_spectra(k, Praw, 200);
            k = k_smooth;
            P = P_smooth;

            loglog(k,abs(P),'linewidth',1,'color','k');

            if E == 1
                loglog(k,0.1*k.^(-1),'k-');
                loglog(k,0.1*k.^(-5/3),'k--');
                loglog(k,0.1*k.^(-7/3),'k-.');
            end

            ind = (k > k_ISR_start) & (k < k_ISR_end);
            if any(ind)
                loglog(k(ind),P(ind),'r'); hold on;
                m(I).phi_uT.ISR(E,S) = get_spectral_exponent(k,P,ind);
            else
                m(I).phi_uT.ISR(E,S) = NaN;
            end

            ind_AE = (k < k_ISR_start);
            if any(ind_AE)
                m(I).phi_uT.AE(E,S) = get_spectral_exponent(k,P,ind_AE);
            else
                m(I).phi_uT.AE(E,S) = NaN;
            end

        
        

        end
    end
    sgtitle(stability_label{I}); 

    subplot(3,5,1); ylabel('\phi_u'); 
    subplot(3,5,6); ylabel('\phi_T');
    subplot(3,5,11); ylabel('\phi_{uT}'); 
    subplot(3,5,1); title(height_labels{1});
    subplot(3,5,2); title(height_labels{2});
    subplot(3,5,3); title(height_labels{3});
    subplot(3,5,4); title(height_labels{4});
    subplot(3,5,5); title(height_labels{5});
end



%% --- Table of ISR slopes by height and stability ---
spec_names = {'phi_uu', 'phi_TT', 'phi_uT'};
spec_labels = {'F_uu', 'F_TT', 'F_uT'};

for sp = 1:3
    fprintf('\n=== %s ISR slopes (|m|) ===\n', spec_labels{sp})
    fprintf('%-14s', 'Height');
    for I = 1:3, fprintf('%18s', stability_label{I}); end
    fprintf('\n%s\n', repmat('-', 1, 68))

    for S = 1:5
        fprintf('z = %.4f m  ', z_vec(S));
        for I = 1:3
            vals = abs(m(I).(spec_names{sp}).ISR(:, S));
            vals = vals(~isnan(vals));
            if numel(vals) >= 2
                fprintf('  %5.2f +/- %4.2f (%d)', mean(vals), std(vals), numel(vals));
            elseif numel(vals) == 1
                fprintf('  %5.2f          (%d)', vals, numel(vals));
            else
                fprintf('       NaN       (0)');
            end
        end
        fprintf('\n');
    end

    % summary excluding lowest height
    fprintf('%-14s', 'z > z_min');
    for I = 1:3
        vals = abs(m(I).(spec_names{sp}).ISR(:, 2:end));
        vals = vals(~isnan(vals));
        fprintf('  %5.2f +/- %4.2f (%d)', mean(vals), std(vals), numel(vals));
    end
    fprintf('\n')
end
