function fig = mtt_fit_eps_Ruddicketal2000(fit_data,varargin)
% INPUT:
%       fit_data: Data structure derived from mtt_fit_eps_Ruddicketal2000.m
%
% optional
%       'figure',nfig : figure handle for data plotting
%       
% OUTPUT:
%       fig: Figure handle
%
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
    if(strcmpi(varargin{i},'verbosity'))
        verbosity = varargin{i + 1};
    end
end


flag_raw = 0;

for i=1:length(varargin)
    if(strcmpi(varargin{i},'figure'))
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
    elseif(strcmpi(varargin{i},'raw'))
         if((i+1) <= length(varargin))
            if(isnumeric(varargin{i+1}))
                k_raw = varargin{i+1}(:,1);
                P_raw = varargin{i+1}(:,2);
                flag_raw = 1;
            end
         end
        
    end
end


legpl = [];
legstr = {};
if(nfig)
    fig = figure(nfig);
else
    fig = figure;
end

% Plot the spectrum
subplot(3,3,[1 2 3 4 5 6])
hold all
% raw data spectrum
if(flag_raw)
    pl = plot(k_raw,P_raw,'-b');
    legpl(end+1) = pl;
    legstr{end+1} = ['raw' ];
end
% original spectrum
pl = plot(fit_data.k,fit_data.P,'-k');
legpl(end+1) = pl;
chilog = floor(log10(fit_data.chi_bestfit));
chistr = sprintf('%.2f',(fit_data.chi_bestfit/10^chilog));
epslog = floor(log10(fit_data.eps_bestfit));
epsstr = sprintf('%.2f',(fit_data.eps_bestfit/10^epslog));
chistrf = ['\chi = ' chistr 'x10^{' num2str(chilog) '}'];
epsstrf = ['\epsilon =' epsstr 'x10^{' num2str(epslog) '}'];
legstr{end+1} = ['int. data \chi' ];

% batchelor fit
pl = plot(fit_data.k,fit_data.P_bat_bestfit,'-r');
legpl(end+1) = pl;
legstr{end+1} = ['fit ( ' chistrf ' ' epsstrf ' )'];
% power law fit
pl = plot(fit_data.k,fit_data.P_pow_bestfit,'-g');
legpl(end+1) = pl;
legstr{end+1} = ['power law fit'];

if(~all(isnan(fit_data.noise)))
    pl = plot(fit_data.k,fit_data.noise);
    legpl(end+1) = pl;
    legstr{end+1} = 'noise';
end
set(gca,'Xscale','log')
set(gca,'Yscale','log')

legend(legpl,legstr,'Location','NorthOutside','Orientation','horizontal')
xlabel('k [cpm]')
ylabel('PSD')

%xlim([floor(min(fit_data.k)),ceil(max(fit_data.k))])

% plot the likelohood matrix
subplot(3,3,7)
hold all
pcolor(log10(fit_data.eps_fit),log10(fit_data.chi_fit),fit_data.C11')
xlabel('log10 \epsilon')
ylabel('log10 \chi')
colorbar
plot(log10(fit_data.eps_bestfit),log10(fit_data.chi_bestfit),'ok')
xlim([min(log10(fit_data.eps_fit(:))),max(log10(fit_data.eps_fit(:)))])
ylim([min(log10(fit_data.chi_fit(:))),max(log10(fit_data.chi_fit(:)))])
%
subplot(3,3,8)
set(gca,'visible','off')
text(.1,.5,['LLR:' num2str(fit_data.LLR)])