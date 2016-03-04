function [chi] = mtt_int_chi(Pxx,k,varargin)
%
%
% Integrates a temperature gradient profile to estimate
% the variance decay of the temperature variance (chi)
%
% input: dTdz, temperature gradient profile [K m^{-1}]
%        fs_k, sampling wavenumber [m^{-1}]
%
% output chi: Temperature variance decay [K^2 s^{-1}]
%
%
% Part of the marine turbulence toolbox:
% https://github.com/MarineDataTools/marine_turbulence_toolbox 

  global mtt_verbosity
  if(~isempty(mtt_verbosity))
    verbosity = mtt_verbosity;
  else
    verbosity = 0;
  end

  for i=1:length(varargin)
    if(strcmpi(lower(varargin{i}),'verbosity'))
      verbosity = varargin{i + 1};
    end
  end

  if(verbosity == 3)
    mtt_message(' ',1)
  end

  dk = k(2) - k(1);
  % The part of the spectrum to integrate and to fit the turbulence spectrum
  k_fit = k;
  P_fit = Pxx;
  Pxx_var = mtt_nansum(Pxx) * dk;
  chi = 6 * 1e-7 * Pxx_var;


