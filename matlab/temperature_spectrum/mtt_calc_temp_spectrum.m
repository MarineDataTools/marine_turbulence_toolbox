function [chi,spec_data] = mtt_calc_temp_spectrum(dTdz,fs_k,varargin)
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

global tt_verbosity
if(~isempty(tt_verbosity))
    verbosity = tt_verbosity;
else
    verbosity = 0;
end

flag_noise = 0;
flag_plot = 0;

for i=1:length(varargin)
    if(strcmpi(varargin{i},'NOISE'))
            size(varargin{i+1})
            noise_k = varargin{i+1}(:,1);
            noise_PSD = varargin{i+1}(:,2);
            flag_noise = 1;
        if(verbosity)            
            tt_message('Found noise data',1)
        end
    else if(strcmpi(varargin{i},'PLOT'))
            flag_plot = 1;
            if((i+1) <= length(varargin))
                if(isnumeric(varargin{i+1}))
                    nfig = varargin{i+1};
                end
            else
                nfig = 0;
            end
            if(verbosity)
                tt_message('Plotting spectrum',1)
            end
        end
    end
end



% The interval in which chi is integrated
k_chi   =  [0.1, 200]; 



% Calculate the spectrum
x = dTdz;
NFFT = length(x);
[Pxx,k] = periodogram(x-nanmean(x),hanning(length(x)),NFFT,fs_k);

% The index of chi integration
ind_chi =  ( k > k_chi(1) ) & ( k < k_chi(2) );

if(flag_noise)
    % Interpolate the noise function onto k
    noise = 10.^(interp1(log10(noise_k),log10(noise_PSD),log10(k),'linear'));
    ind_noise = find(Pxx < noise);
    if(~isempty(ind_noise))
        ind_noise = min(ind_noise);
        % Now we want a half a decade of more noise
        k_noise = k(ind_noise);
        ind_some_noise = find(log10(k) > (log10(k_noise) + 0.5));
        if(~isempty(ind_some_noise) > 0)
            ind_some_noise = min(ind_some_noise);
        else
            ind_some_noise = ind_noise;
        end
    end
    
    ind_chi(ind_some_noise:end) = 0;
    
    Pxx_denoise = Pxx - noise;
else
    Pxx_denoise = Pxx;
end



dk = k(2) - k(1);
% The part of the spectrum to integrate and to fit the turbulence spectrum
k_fit = k(ind_chi);
P_fit = Pxx(ind_chi);
Pxx_var = nansum(Pxx_denoise(ind_chi)) * dk;
chi = 6 * 1e-7 * Pxx_var;

spec_data.k = k;
spec_data.P = Pxx;
spec_data.ind_chi = ind_chi;
spec_data.chi = chi;
if(flag_noise)
    spec_data.noise = noise;
end


% ------------------------ Plot the data ----------------------------------
if(flag_plot)
    legpl = [];
    legstr = {};
    if(nfig)
        figure(nfig)
    else
        figure
    end
    hold all
    
    pl = plot(k,Pxx,'-k');
    legpl(end+1) = pl;
    legstr{end+1} = 'orig. data';
    pl = plot(k(ind_chi),Pxx(ind_chi),'-r');
    legpl(end+1) = pl;
    
    chilog = floor(log10(chi));
    chistr = sprintf('%.2f',(chi/10^chilog));
    
    legstr{end+1} = ['int. data ( \chi=' chistr 'x10^{' num2str(chilog) '} )'];
    if(flag_noise)
        %plot(k,noise)
        pl = plot(noise_k,noise_PSD);
        legpl(end+1) = pl;
        legstr{end+1} = 'noise';
    end
    set(gca,'Xscale','log')
    set(gca,'Yscale','log')
    legend(legpl,legstr,'Location','NorthOutside','Orientation','horizontal')
    xlabel('k [cpm]')
    ylabel('PSD')
    xlim([floor(min(k)),ceil(max(k))])
    
end
    


return
