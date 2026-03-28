function fit = fit_F1_coeff(k, Pk, coeff_set)
% for large scales 
k = k(:); Pk = Pk(:);
mdl = @(q,k) q(1).*coeff_set.*k.^(-1);

% Initial guesses and bounds
C = 0.5;      
lb = 0;    
ub = Inf;

% Log-space residuals for scale invariance
obj  = @(q) log10(abs(mdl(q,k))) - log10(abs(Pk));
opts = optimoptions('lsqnonlin','Display','off','MaxIter',3000, ...
                    'FunctionTolerance',1e-12,'StepTolerance',1e-12);
q = lsqnonlin(obj,C,lb,ub,opts);

C = q(1); 
Phat = mdl(q,k);
fit = struct('C',C,'k',k,'Phat',Phat);
% quick plot
%figure;
%loglog(kappa, abs(y), '.', kappa, abs(yhat), '-', 'LineWidth',1.8);
%xlabel('\kappa = k z'); ylabel( '|k B(k)|');
%grid on; legend('data','Kaimal fit');
%title(sprintf('A=%.3g, D=%.3g, p=%.3g  (high-k slope → -%.3g)',A,D,pexp,pexp));
end
