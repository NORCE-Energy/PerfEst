function [tapering_coefficient,output] = tapering_function_theoretical(corr_field,ensize,varargin)
% function tapering_coefficient = tapering_function(corr_field,,corr_field_shuffled,varargin)
%
% -- INPUTS --
% corr_field:               correlation fields
% ensize:                   ensemble size
% varargin:                 reserved for future development
%
% -- OUTPUTS --
% tapering_coefficient:     tapering coefficients
% output:                   other information
%
% Made by XILU, Jan 2018
% Copyright (c) 2010-2018 IRIS, All Rights Reserved.
% $Id:  $
% $DateTime:  $

output = [];

fh = @(x,N) (N-1).* x.^2 ./ (x.^4 + (N-3) .* x.^2 + 1);

tapering_coefficient = fh(corr_field,ensize);
