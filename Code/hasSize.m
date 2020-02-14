function r = hasSize(array, dim)
% function r = hasSize(array, dim)
%
% Helper function to check if array has given size. It handles the
% case where sizes are compatible but not identical like
% ```
%     [ 2, 3, 1 ] == [ 2, 3 ]
% ```
%
% Inputs:
% -------
% array:       The array to check the size of.
% dim:         The size to check for.
%
% Outputs:
% --------
% r:           True if the size is compatible.
%
% Copyright (c) 2016-2017 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/anma/mlabScripts/hasSize.m#1 $
% $DateTime: 2017/12/19 16:10:07 $

a_dim = size(array);
a_ndims = ndims(array);

r = isequal(a_dim, dim(1:a_ndims)) && all(dim(a_ndims + 1:end) == 1);

end
