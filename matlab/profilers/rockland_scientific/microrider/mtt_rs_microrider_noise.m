function noise = tt_rockland_microrider_noise(k)
%
%
% noise spectrum of the rockland microrider temperature amplifier
% circuits, functional form taken from
% Luc Rainville and Peter Winsor: "Mixing across the Arctic Ocean:
% Microstructure observations during the Beringia 2005 Expedition",
% GRL 2008
%    
% input k: Wavenumber [cpm]
% output noise: spectral power density [(K m^{-1})^2/cpm]
%   
%
% part of the turbulence toolbox
    
    EL_NOISE_k  =  logspace(0,5.2);
    EL_NOISE    =  logspace(log10(2e-9),-2.69);  
    
    noise = interp1(log10(EL_NOISE_k),log10(EL_NOISE),log10(k));
    noise = 10.^noise;
    
