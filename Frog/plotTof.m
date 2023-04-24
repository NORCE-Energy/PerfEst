function h = plotTof

load trueSolutionSmoother.mat;
load inputData.mat;
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
h = figure('Name','Frog time of flight','NumberTitle','off',...
       'position',[1, 776, 1165, 335]);
subplot(1,3,1);
v = tofArtObserved.*actnumMask;
b = imagesc(v);
set(b,'AlphaData',~isnan(v))
title('TOF Arterial True');
colorbar;
cl = get(gca,'clim');
set(gca,'fontsize',14,'fontweight','bold')
subplot(1,3,2);
v = state.tofArt.*actnumMask;
b=imagesc(v);
set(gca,'clim',cl);
set(b,'AlphaData',~isnan(v))
title('TOF Arterial Final');
colorbar;
set(gca,'fontsize',14,'fontweight','bold')
subplot(1,3,3);
v = (state.tofArt-tofArtObserved).*actnumMask;
b = imagesc(v);
set(gca,'clim',cl);
set(b,'AlphaData',~isnan(v))
title('TOF Arterial Diff')
colorbar;
set(gca,'fontsize',14,'fontweight','bold')
drawnow;