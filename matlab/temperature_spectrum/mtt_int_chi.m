function [chi] = mtt_int_chi(Pxx,k,varargin)
%
% Integrates a temperature gradient profile to estimate
% the variance decay of the temperature variance (chi)
% input: dTdz, temperature gradient profile [K m^{-1}]
%        fs_k, sampling wavenumber [m^{-1}]
%        noise (Optional), vector of noise function
%               e.g. tt_int_chi(dTdz,fs_k,'noise',noise)
%               noise must have the size (n,2) where n is the number of
%               data records; this is an example of electronic noise record
%               for a Rockland microrider:
%               EL_NOISE_k  =  logspace(0,5.2);
%               EL_NOISE    =  logspace(log10(2e-9),-2.69);
%               noise = [EL_NOISE_k',EL_NOISE_k'];
%        plot  (Optional), plot the spectrum, the integration limits and
%               other useful information
%               e.g. tt_int_chi(dTdz,fs_k,'noise',noise,'plot')
% output chi: Temperature variance decay [K^2 s^{-1}]
%
%
% this is a part of the turbulence toolbox
% 

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


