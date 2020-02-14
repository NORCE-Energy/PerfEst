function H = plotPerf(varargin)


fs = setProperty(varargin,'fontSize',14);
fw = setProperty(varargin,'fontWeight','bold');
layer = setProperty(varargin,'layer',1);

load trueSolutionSmoother.mat;
load inputData.mat;
if exist('initialState.mat','file')
    load('initialState.mat','options','initialState','ref','prm');
    initialOptions = options; %#ok<*NASGU>
end
if exist('finalState.mat','file')
    load finalState.mat;
end
if ~exist('state','var') 
    if exist('initialState','var')
        state = initialState;
    else
        error('No results to plot.');
    end
end

actnumMask = reshape(options.actnum,options.dim(1),options.dim(2),options.dim(3));
actnumMask(actnumMask==0) = NaN;
H = figure('Name','Tissue perfusion','NumberTitle','off',...
       'position',[1, 776, 720, 284]);
h(1) = subplot(1,3,1);
Pcrs(initialState.rateQ<1e-6)=0;
uc = 6000; % unit conversion mL/min/100mL
v = Pcrs(:,:,layer).*actnumMask(:,:,layer)*uc;
b = imagesc(v);
set(b,'AlphaData',~isnan(v))
title('True','fontSize',fs,'fontweight',fw);
%colorbar;
set(gca,'xtick',[])
set(gca,'ytick',[])
%set(gca,'fontsize',fs,'fontweight',fw)
set(h(1),'units','centimeters','position',[1.5 1.5 4 4.8484]);
h(2) = subplot(1,3,2);
perf_ini = initialState.rateQ;
perf_ini(initialState.rateQ<1e-6)=0;
v = perf_ini(:,:,layer).*actnumMask(:,:,layer)*uc;
b=imagesc(v);
set(gca,'xtick',[])
set(gca,'ytick',[])
set(b,'AlphaData',~isnan(v))
cl = get(gca,'clim');
title('Prior','fontSize',fs,'fontweight',fw);
%colorbar;
%set(gca,'fontsize',fs,'fontweight',fw)
set(h(2),'units','centimeters','position',[6.5 1.5 4 4.8484]);
h(3) = subplot(1,3,3);
%v = (Pcrs-perf_end).*actnumMask;
perf_end = state.rateQ;
perf_end(initialState.rateQ<1e-6)=0;
v = perf_end(:,:,layer).*actnumMask(:,:,layer)*uc;
b = imagesc(v);
set(b,'AlphaData',~isnan(v))
title('Posterior','fontSize',fs,'fontweight',fw)
colorbar;
%set(gca,'fontsize',fs,'fontweight',fw)
set(h(3),'units','centimeters','position',[11.5 1.5 4 4.8484]);
set(gca,'xtick',[])
set(gca,'ytick',[])

for I = 1:3
    set(h(I),'clim',cl);
end
