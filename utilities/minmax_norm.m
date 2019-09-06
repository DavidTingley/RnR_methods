function [m] = minmax_norm(mm)

m = mm - min(mm);
m = m ./ max(m);
m(isnan(m)) = 0;
return 