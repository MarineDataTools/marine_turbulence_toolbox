function [Pxx ,k , Pxx_denoise, Pxx_noise] = mtt_calc_spectrum(x,fs,varargin)
%
% 
% input: 
%        fs_k, sampling wavenumber [m^{-1}]
%        noise (Optional), vector of noise function
%               e.g. tt_int_chi(dTdz,fs_k,'noise',noise)
%               noise must have the size (n,2) where n is the number of
%               data records; this is an example of electronic noise record
%               for a Rockland microrider:
%               EL_NOISE_k  =  logspace(0,5.2);
%               EL_NOISE    =  logspace(log10(2e-9),-2.69);
%               noise = [EL_NOISE_k',EL_NOISE_k'];
% output:
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

if(verbosity == 3)
  mtt_message(' ',1);
end

flag_noise = 0;
flag_plot = 0;

% Read in local verbosity
for i=1:length(varargin)
  if(strcmpi(varargin{i},'verbosity'))
    verbosity = varargin{i+1};
  end
end

% Read in local verbosity
flag_hanning = 0;
for i=1:length(varargin)
  if(strcmpi(varargin{i},'hanning'))
    flag_hanning = 1;
  end
end

for i=1:length(varargin)
    if(strcmpi(varargin{i},'NOISE'))
            noise_k = varargin{i+1}(:,1);
            noise_PSD = varargin{i+1}(:,2);
            flag_noise = 1;
        if(verbosity)            
            mtt_message('Found noise data',1)
        end
    end
end


NFFT = length(x);
if(flag_hanning)
    NFFT = floor(length(x));
    [Pxx,k] = periodogram(x-mtt_nanmean(x),hanning(length(x)),NFFT,fs);
    Xd = detrend(x,'linear');
    %figure(100)
    %hold all
    %plot(x)
    %plot(Xd)
    %pause
    %[Pxx, k] = pwelch(x-mtt_nanmean(x), [], [], [], fs); % Raw Power Spectrum
    %[Pxx, k] = pwelch(Xd, [], [], [], fs); % Raw Power Spectrum
    dk = k(2) - k(1);
    %Pxx = Pxx / sum(Pxx*dk) * var(x);    
    var1 = var(x);
    var2 = sum(Pxx*dk);
    var3 = 1/length(x) * sum((x - mean(x)).^2);

    disp([' calc spectrum var1: ' num2str(var1) ' var2: ' num2str(var2) ' var3: ' num2str(var3)])
    %[Pxx, k] = pwelch(x, [], [], [], fs); % Raw Power Spectrum
    Pxx = Pxx / sum(Pxx*dk) * var(x);
else
    %[Pxx,k] = periodogram(x-mtt_nanmean(x),[],NFFT,fs);
end


if(flag_noise)
    % Interpolate the noise function onto k
    noise = 10.^(interp1(log10(noise_k),log10(noise_PSD),log10(k),'linear'));
    Pxx_denoise = Pxx - noise;
    Pxx_noise = noise;  
else
    Pxx_denoise = Pxx * NaN;
    Pxx_noise = Pxx * NaN;  
end



