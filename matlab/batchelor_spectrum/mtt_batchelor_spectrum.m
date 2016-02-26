function P = mtt_batchelor_spectrum(k,Eps,Chi,nu,nut)
% Computes the one-dimensional Batchelor spectrum for small-scale 
% temperature gradients as function of the cyclic wavenumber k. Two
% different algorithms exists based on the direct evaluation of the
% integrals and a polynomial expansion (see Luketina and Imberger 2001),
% respectively. 
%
% INPUT PARAMTERS:
% k:                    wavenumber vector in cyclic units (cpm)
% Eps:                  dissipation rates to fit spectrum against (W/kg)
% Chi:                  mixing rate of temperature variance to fit spectrum against (degC^2/s) 
% nu:                   viscosity (m^2/s)
% nut:                  diffusivity of temperature (m^2/s)
% 
% OUTPUT PARAMETERS:    
% P:                    vector of Batchelor spectrum (degC^2/m^2*cpm)
%
% 
% Part of the marine turbulence toolbox:
% https://github.com/MarineDataTools/marine_turbulence_toolbox


% define model constants
q  = 3.4;
c  = sqrt(2*q);

% Batchelor wavenumber (cyclic units)
kb = 1/2/pi*(Eps/nu/nut/nut)^0.25;

% non-dimensional wavenumber
a  = c*k/kb;

% compute Batchelor spectrum
% F  = a .* ( gauss(a) - a .* direct_int(a) );
F  = a .* ( gauss(a) - a .* expand_int(a) );
P  = 0.5*c*Chi/kb/nut*F;

% check for negative values (numerical artefacts)
PMin     = min(P(P>0));
if(isempty(PMin)) % No data at all
    P(P<=0)  = 0;
else
    P(P<=0)  = PMin;
end

end

function y = gauss(x)
% Simple Gauss kernel

y = exp(-0.5 * x .* x);

end

function y = direct_int(x)
% Directly compute integral appearing in Eq. (9) of Luketina and Imberger
% (2001). The algorithm used here is based on the following decomposition:
% Int_x^inf = Int_0^inf - Int_0^x1 - Int_x1^x

Int0 = sqrt(2*pi)/2;           % Int_0^infty (exp(-x^2/2) dx
Int1 = quad(@gauss,0,x(1));
y    = Int0 - Int1 -  cumtrapz(x,gauss(x));

end

function y = expand_int(x)
% Alternatively, use the integral expansion suggested in Eqs. (12) and (13) 
% by Luketina and Imberger (2001). This expansion is very accurate and 
% usually faster.
 
% model constants
a1 = 0.319381530;
a2 = 0.356563782;
a3 = 1.781477937;
a4 = 1.821255978;
a5 = 1.330274429;
a6 = 0.2316419;

t = 1 ./ (1 + a6*x);
y = gauss(x) .* (a1*t - a2*t.*t + a3*t.*t.*t - a4*t.*t.*t.*t + a5*t.*t.*t.*t.*t);

end




