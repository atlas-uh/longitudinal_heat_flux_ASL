function [slope] = getm(k,P,ind)


x = log10(k(ind));
y = log10(P(ind));
is_inf_y = isinf(y);
finite_indices = find(~is_inf_y);
p = polyfit(x(finite_indices), y(finite_indices), 1);    % linear fit
slope = p(1);          % slope in log-log
            


