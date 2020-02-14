function H=eplot(M,G)
% eplot -  Plot matrix values centered on face, not on grid
%
% function H=eplot(M,G)
%
% Designed for eclipse numbering (using z-axis), i.e. 
% starting with first element at top layer.
% (eclipse numbering runs first in x, then y, then z direction).
% 
% The function is based on pplot.
%
% This function makes a 2D plot of the elements of M, using the
% (normalized) value of each element as an index into the color
% map.  The plot is face centered, i.e., each matrix element
% corresponds to a grid cell, and not to a vertix.
% The function defaults to plotting the grid lines, but these may
% be turned off (better viewing for large matrices) by setting
% G = 0.

% Copyright (c) 2006 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/main/plot2d/eplot.m#2 $
% $DateTime: 2014/09/19 17:08:28 $

%
if nargin==1
  T=pplot(M');
elseif nargin==2
  T=pplot(M',G);
else
  error('Number of argument to function eplot should be 1 or 2')
end
set(gca,'ydir','reverse');
axis image;
if nargout==1
  H=T;
end
