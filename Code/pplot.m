function H=pplot(M,G)

%pplot -  Plot matrix values centered on face, not on grid
%
% function H=pplot(M,G)
%
% This function makes a 2D plot of the elements of M, using the
% (normalized) value of each element as an index into the color
% map.  The plot is face centered, i.e., each matrix element
% corresponds to a grid cell, and not to a vertix.
% The function defaults to plotting the grid lines, but these may
% be turned off (better viewing for large matrices) by setting
% G = 0.

% Copyright (c) 2006 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/main/plot2d/pplot.m#1 $
% $DateTime: 2011/10/13 14:36:15 $

%
newplot ;
[a,b]=size(M);
zfill = mean(M(1:a,b));
x = (1:b+1)-.5;
y = (1:a+1)-.5;
%
N = [ M          , ones(a,1)*zfill;...
      ones(1,b)*zfill , zfill ] ;
%
S=pcolor(N);
set(S,'xdata',x,'ydata',y);
set(gca,'xlim',[.5,b+.5],'ylim',[.5,a+.5]);
%
if nargin > 1 && G == 0
  shading('flat');
end
if nargout == 1
  H=S;
end
%
