function [y,ind_max] = mtt_nanmax(x)


  ind = find(~isnan(x));
  [y,ind_max_tmp] = max(x(ind));
  ind_max = ind(ind_max_tmp);
