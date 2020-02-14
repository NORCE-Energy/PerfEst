function property = setProperty(v,propStr,default)

% function property = setProperty(v,propStr,default)
%
% set default values to variable input arguments (varargin)
% --
% v:        varargin
% propStr:  name(s) of the variable input arguments
% default:  default values
% Example usage: timeShift = setProperty(varargin,'timeShift',0) in "plotForecast.m"
%
% Copyright(c) International Research Institute of Stavanger (IRIS)
% $Id: //depot/rfmatlab/main/iofun/setProperty.m#3 $
% $DateTime: 2017/09/29 11:25:04 $


vn = v(1:2:end);
ind = cellfun(@(x) ~isempty(strfind(x,propStr)),vn);
ind = find(ind == 1,1);
if ~isempty(ind)
    property = v{2*ind};
else
    property = default;
end