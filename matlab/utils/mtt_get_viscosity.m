function vis = mtt_get_viscosity(T)
%
% mtt_get_viscosity(T)
% input:
% T [degC]
% output:
% viscosity [m**2/s]
% vis = (1.792747-(.05126103*T)+(0.0005918645*T*T))*1e-6;
%
% Peter Holtermann adapted from Ilker Fer
%
% Part of the marine turbulence toolbox:
% https://github.com/MarineDataTools/marine_turbulence_toolbox    

vis = (1.792747 - (.05126103 * T) + (0.0005918645 * T * T )) * 1e-6;
