function fit = fit_Production(k, Pk)

k = k(:); Pk = Pk(:);

% Model (A,D,p free); bake in the sign so log-residuals are well-defined
mdl = @(q,x) (q(1)./x) .* (1 + q(2).*x.^2).^-q(3); % kB model

% Initial guesses and bounds
A0 = median(abs(Pk));
B0 = 30;                 % Kaimal knee O(10–40)
beta0 = 0.5;                % start near 5/3; free to move
lb = [0     0   0.1];    % A>=0 (sign handled by 'sgn'), D>=0, p>=0.5
ub = [Inf   Inf 10];

% Log-space residuals for scale invariance
obj  = @(q) log10(abs(mdl(q,k))) - log10(abs(Pk));
opts = optimoptions('lsqnonlin','Display','off','MaxIter',3000, ...
                    'FunctionTolerance',1e-12,'StepTolerance',1e-12);
q = lsqnonlin(obj,[A0 B0 beta0],lb,ub,opts);

A = q(1); B = q(2); gamma = q(3);
Pkhat = mdl(q,k);

fit = struct('A',A,'B',B,'gamma',gamma,'k',k,'Pk',Pk,'Pkhat',Pkhat);
% quick plot
%figure;
%loglog(kappa, abs(y), '.', kappa, abs(yhat), '-', 'LineWidth',1.8);
%xlabel('\kappa = k z'); ylabel( '|k B(k)|');
%grid on; legend('data','Kaimal fit');
%title(sprintf('A=%.3g, D=%.3g, p=%.3g  (high-k slope → -%.3g)',A,D,pexp,pexp));
end
