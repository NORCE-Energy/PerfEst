function W = concatenateW(Wa, Wb)
% function W = concatenateW(Wa, Wb)
%
% Concatenate matrix as either vector or block diagonal matrix.
%
% Vectors represents diagonal matrices with the given elements
% on the diagonal. If both inputs are vectors, they are simply
% concatenated, otherwise they are converted to a block diagonal
% matrix.
%
% If an input is diagnoal (as indicated by "isdiag" function), it
% is treated as a vector.
%
% INPUTS
% ------
% Wa:        First matrix in concatination
% Wb:        Second matrix in concatination
%
% OUTPUTS
% -------
% W:         Concatination of the two inputs
%
% Copyright (c) 2017 IRIS, All Rights Reserved.
% $Id: $
% $DateTime: $

if (size(Wa, 2) == 1 || isdiag(Wa)) && (size(Wb, 2) == 1 || isdiag(Wb))
    % Concatenate as vectors

    if size(Wa, 2) ~= 1
        Wa = diag(Wa);
    end

    if size(Wb, 2) ~= 1
        Wb = diag(Wb);
    end

    W = [ Wa; Wb ];
else
    % Concatenate as blocks

    if size(Wa, 2) == 1
        Wa = diag(Wa);
    end

    if size(Wb, 2) == 1
        Wb = diag(Wb);
    end

    W = blkdiag(Wa, Wb);
end

end
