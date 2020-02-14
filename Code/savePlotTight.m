function savePlotTight(prtFile,h,varargin)

% save plots into different formats (default = fig, eps, png)
% 
% Copyright (c) 2012-2016 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/main/graphics/savePlotTight.m#3 $
% $DateTime: 2018/12/14 13:49:52 $

save_fig = setProperty(varargin,'save_fig',1);
save_eps = setProperty(varargin,'save_eps',1);
save_png = setProperty(varargin,'save_png',1);

if save_fig && nargin > 1   
    saveas(h,[prtFile '.fig']);
end

if save_eps
    print('-r0','-depsc2',prtFile);
end

if save_png
    print('-r0','-dpng',prtFile);
    system(['convert ',prtFile,'.png -trim ',prtFile,'.png']);
end
