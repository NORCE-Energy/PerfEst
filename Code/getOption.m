function value = getOption(options,key,default)

% function value = getOption(options,'key',default)
%
% Get option if it exist, or return NaN
%
% Input:
% ------
% options : structure with options
% key     : structure key
% default : return this value if key is not present 
%
% Output: 
% -------
% value   : structure value
%
% Copyright(c) International Research Institute of Stavanger (IRIS)
% $Id: //depot/rfmatlab/main/iofun/getOption.m#1 $
% $DateTime: 2017/09/29 11:24:28 $


value = NaN;
if isfield(options,key) 
    value  = getfield(options,key);
elseif nargin == 3
    value = default;
end