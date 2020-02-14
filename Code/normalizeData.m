function normData = normalizeData(data, sW)
% function normData = normalizeData(data, sW)
%
% Normalizes measurement data using given std. The std
% can be eiter a vector of individual stds, or a matrix
% (typically transposed Cholesky factor of covariance).
%
% INPUTS
% ------
% data:          Data to be normalized
% sW:            Std of the data, as either vector or matrix
%
% OUTPUTS
% -------
% normData:      Normalized data

if size(sW, 2) == 1
    normData = data ./ repmat(sW, 1, size(data, 2));
else
    normData = sW\data;
end

end
