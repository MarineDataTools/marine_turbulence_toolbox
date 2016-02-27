% This script processes data from a Rockland Scientific Microrider
%
% 
% 
% Part of the marine turbulence toolbox:
% https://github.com/MarineDataTools/marine_turbulence_toolbox      

clear
close all

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

global mtt_verbosity
mtt_verbosity = 3;

% This is for octave (TODO version check)
if(isOctave)
    struct_levels_to_print(0)
end

% Loading a temperature microstructure profile
load(['DAT_310p.mat']);
who

% We need temperature, pressure and time
disp 'Loading data'
T = uMnc.T1;
p = uMnc.p;
t = uMnc.time;

profile.T1 = uMnc.T1;
profile.T2 = uMnc.T2;
profile.p  = uMnc.p;
profile.time  = uMnc.time;
profile.lon  = 15.9882;
profile.lat  = 55.2518;
profile.name = 'DAT_310p.mat';
profile.date = datenum(2009,1,1);

mtt_plot_profile(profile,'verbosity',1);
disp 'Finished'

disp 'Calculating velocity'
[w] = mtt_calc_gradient(p,t);
disp 'Calculate temperature gradient'
[dTdp] = mtt_calc_gradient(T,p);


% Noise function of the microrider
EL_NOISE_k  =  logspace(0,5.2);
EL_NOISE    =  logspace(log10(2e-9),-2.69);

eps_fit = logspace(-12,-5,50);
chi_fit = logspace(-12,-5,50);
% Calculate the turbulence in depths intervals of dp, between
% max(p_turb) and min(p_turb)
dp = 1.0;
p_turb12 = 25:dp:70;
p_turb = p_turb12(1:end-1) + dp/2;
for i=1:length(p_turb)
    ind_data = (p > p_turb12(i)) & (p <= p_turb12(i+1));
    if(sum(ind_data) > 100)
        disp(['Data in ' num2str(p_turb(i))])
        figure(2)
        clf
        subplot(1,2,1)
        plot(T(ind_data),p(ind_data))
        axis ij
        subplot(1,2,2)
        plot(w(ind_data),p(ind_data))
        axis ij
        
        t_seg = t(ind_data);
        T_seg = T(ind_data);  
        dTdp_seg = dTdp(ind_data);
        w_seg = w(ind_data);
        w_seg_avg = mtt_nanmean(w_seg);
        fs_t = 1/(t_seg(2) - t_seg(1));
        fs_k = abs(fs_t/w_seg_avg);
        ind_nan = ~isnan(dTdp_seg);
        if(sum(ind_nan) > 25)
            figure(3)
            clf
            figure(4)
            clf
            [Pxx,k,Pxx_denoise,Pxx_noise] = mtt_calc_spectrum(dTdp_seg(ind_nan),fs_k,'noise',[EL_NOISE_k',EL_NOISE']);
            % The interval in which chi is integrated
            k_chi   =  [0.1, 200];
            % The index of chi integration
            ind_chi =  ( k > k_chi(1) ) & ( k < k_chi(2) );
            [chi(i)] = mtt_int_chi(Pxx_denoise(ind_chi),k(ind_chi));
            vis = mtt_get_viscosity(mtt_nanmean(T(ind_data)));
            mtt_fit_eps_Ruddicketal2000(k(ind_chi),Pxx(ind_chi),chi_fit,eps_fit,vis,4,'noise',Pxx_noise(ind_chi),'plot',4)
            pause
        else
            disp('not enough data')
        end
        
    end
end


