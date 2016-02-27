function y = mtt_nansum(x)


  ind = ~isnan(x);
  y = sum(x(ind));
