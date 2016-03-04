function y = mtt_nanmean(x)


  ind = ~isnan(x);
  y = mean(x(ind));
