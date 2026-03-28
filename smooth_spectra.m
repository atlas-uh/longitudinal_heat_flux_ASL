function [f,output,slope]=smooth_spectra(F,input,bins)
A=[F input];
output2=zeros(bins,size(A,2));
slope = zeros(bins,1);
binedge=logspace(min(log10(A(2:end,1))),max(log10(A(2:end,1))),bins+1);
warning off
for ii=1:bins
    binlow=binedge(ii);
    binhigh=binedge(ii+1);
    temp=A(((A(:,1)>=binlow)&(A(:,1)<binhigh)),:);
    output2(ii,:)=mean(temp,1);
   
    x = log10(temp(:,1));
    y = log10(temp(:,2));
    is_inf_y = isinf(y);
finite_indices = find(~is_inf_y);
x = x(finite_indices);
y = y(finite_indices);
is_inf_x = isinf(x);
finite_indices = find(~is_inf_x);
x = x(finite_indices);
y = y(finite_indices);
p = polyfit(x, y, 1);    % linear fit

   
    slope(ii) = p(1);          % slope in log-log
end
warning on
output=output2(isfinite(output2(:,1)),:);
slope=slope(isfinite(output2(:,1)));
f=output(:,1);
output=output(:,2);