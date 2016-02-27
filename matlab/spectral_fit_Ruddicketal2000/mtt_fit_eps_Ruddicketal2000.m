function [chi_bestfit,eps_bestfit,LLR] = mtt_fit_eps_Ruddicketal2000(k,P,chi_fit,eps_fit,vis,d,varargin)
%
% Fits the dissipation rate epsilon 
% using the algorithm published in 
% "Maximum Likelihood Spectral Fitting: The Batchelor Spectrum", Ruddick
% B., Anis A., Thompson K., JOURNAL OF ATMOSPHERIC AND OCEANIC TECHNOLOGY,
% Vol. 17, pg. 1541-1555, 2000.
%
% INPUT:
%       k: wavenumber [cycles per meter]
%       P: Power spectrum [K^2 m^{-1}]
%       chi_fit: Temperature variance decay rate to fit spectrum against
%       eps_fit: Dissipation rate of turbulent kinetic energy  to fit spectrum against [W kg^{-1}]
%       vis: Viscosity [m^2 s^{-1}]
%       d: Degree of freedom of power spectrum
% optional
%       'noise',noise:
%       'plot',[fig]
%
% OUTPUT:
%
% Part of the marine turbulence toolbox:
% https://github.com/MarineDataTools/marine_turbulence_toolbox      

global mtt_verbosity
if(~isempty(mtt_verbosity))
    verbosity = mtt_verbosity;
else
    verbosity = 0;
end

flag_plot = 0;
flag_noise = 0;

for i=1:length(varargin)
    if(strcmpi(varargin{i},'NOISE'))
        if(verbosity)
            mtt_message('Found noise data',1)
        end
            if((size(varargin{i+1},1) == length(k)) && (size(varargin{i+1},2) == 1))
                if(verbosity>1)
                    mtt_message('Noise data has the same length as k and only one dimension, taking it without interpolating',1)
                end
                noise = varargin{i+1}(:);
            else
                noise_k = varargin{i+1}(:,1);
                noise_PSD = varargin{i+1}(:,2);
            end
            
            flag_noise = 1;

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
                mtt_message('Plotting spectrum',1)
            end
        end
    end
end


dk = k(2) - k(1);


% The epsilon to fit against
nut = 1e-7; % The diffusivity of temperature

kb_bat = ((eps_fit/vis/vis/nut).^0.25)/(2*pi); % Batchelor wavenumber
P_fit_bat = zeros(length(k),length(eps_fit),length(chi_fit)); % preallocate array for the batchelor spectra
P_th = zeros(size(P_fit_bat)); % preallocate array for the theoretical spectra (Batchelor plus noise)
C11 = zeros(length(eps_fit),length(chi_fit)); % preallocate array for the likelihood fit
for l=1:length(eps_fit)
    for m=1:length(chi_fit)
        P_fit_bat(:,l,m) = mtt_batchelor_spectrum(k,eps_fit(l),chi_fit(m),vis,nut);
        if(flag_noise)
            P_th(:,l,m) = P_fit_bat(:,l,m) + noise;
        else
            P_th(:,l,m) = P_fit_bat(:,l,m);
        end
        pdffac1 = ( d * P ) ./ ( P_th(:,l,m) );
        C11i = d./(P_th(:,l,m)) .* chi2pdf(pdffac1 , d);
        logC11i = log(C11i);
        logC11i(isinf(logC11i)) = NaN;
        C11(l,m)  = mtt_nansum(logC11i);
    end
end

C11(isinf(C11)) = NaN;

% Find the maximum likelihood
[~,ind_max] = mtt_nanmax(exp(C11(:)/mtt_nanmax(C11(:))));
[ind_eps,ind_chi] = ind2sub(size(C11),ind_max);
eps_bestfit = eps_fit(ind_eps);
chi_bestfit = chi_fit(ind_chi);
P_bat_bestfit = mtt_batchelor_spectrum(k,eps_bestfit,chi_bestfit,vis,nut) + noise;
C11_max = C11(ind_max);


% Fit a line in log-log space of the form S1_th = Ak^(-b) + S1_n
b = -1:0.05:10;
for l=1:length(b)
    P_pow = k.^(-b(l)) + noise; % The spectrum
    A = (mtt_nansum(P) * dk) / (mtt_nansum(P_pow) * dk);
    P_pow = A * P_pow;
    pdffac1 = ( d * P ) ./ ( P_pow );
    C11i_pow = d./(P_pow) .* chi2pdf(pdffac1 , d);
    logC11i_pow = log(C11i_pow);
    logC11i_pow(isinf(logC11i_pow)) = NaN;
    C11_pow(l)  = mtt_nansum(logC11i_pow);
end

[~,ind_pow] = max(exp(C11_pow(:)/max(C11_pow)));
P_pow = k.^(-b(ind_pow)) + noise;
A = ( mtt_nansum(P) * dk ) / ( mtt_nansum(P_pow) * dk );
P_pow_bestfit = A * k.^(-b(ind_pow)) + noise;
C11_pow_max = C11_pow(ind_pow);

LLR = log10(exp(C11_max-C11_pow_max));



% ------------------------ Plot the data ----------------------------------
if(flag_plot)
    legpl = [];
    legstr = {};
    if(nfig)
        figure(nfig)
    else
        figure
    end
    
    % Plot the spectrum
    subplot(3,3,[1 2 3 4 5 6])
    hold all
    % original spectrum
    pl = plot(k,P,'-k');
    legpl(end+1) = pl;
    chilog = floor(log10(chi_bestfit));
    chistr = sprintf('%.2f',(chi_bestfit/10^chilog));
    epslog = floor(log10(eps_bestfit));
    epsstr = sprintf('%.2f',(eps_bestfit/10^epslog));
    chistrf = ['\chi = ' chistr 'x10^{' num2str(chilog) '}'];
    epsstrf = ['\epsilon =' epsstr 'x10^{' num2str(epslog) '}'];
    legstr{end+1} = ['int. data \chi' ];
    
    % batchelor fit
    pl = plot(k,P_bat_bestfit,'-r');
    legpl(end+1) = pl;
    legstr{end+1} = ['fit ( ' chistrf ' ' epsstrf ' )'];
    % power law fit
    pl = plot(k,P_pow_bestfit,'-g');
    legstr{end+1} = ['power law fit'];
    
    for l=1:10:length(eps_bestfit)
        plot(k,P_th(:,l),'-','Color',[.5, .5, .5])
    end
    
    if(flag_noise)
        pl = plot(k,noise);
        legpl(end+1) = pl;
        legstr{end+1} = 'noise';
    end
    set(gca,'Xscale','log')
    set(gca,'Yscale','log')
    legend(legpl,legstr,'Location','NorthOutside','Orientation','horizontal')
    xlabel('k [cpm]')
    ylabel('PSD')
    xlim([floor(min(k)),ceil(max(k))])
    
    % plot the likelohood matrix
    subplot(3,3,7)
    hold all
    pcolor(log10(eps_fit),log10(chi_fit),C11')
    xlabel('log10 \epsilon')
    ylabel('log10 \chi')
    colorbar
    plot(log10(eps_fit(ind_eps)),log10(chi_fit(ind_chi)),'ok')
    xlim([min(log10(eps_fit(:))),max(log10(eps_fit(:)))])
    ylim([min(log10(chi_fit(:))),max(log10(chi_fit(:)))])
    % 
    subplot(3,3,8)
    set(gca,'visible','off')
    text(.1,.5,['LLR:' num2str(LLR)])
end
    


return
