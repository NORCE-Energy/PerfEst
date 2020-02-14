function outW = sqrtW(inW)
% function outW = sqrtW(inW)
%
% Returns eitehr the square root or transposed
% cholesky factor of input. Typically used to
% convert variance to std.


if size(inW, 2) == 1
    outW = sqrt(inW);
else
    outW = chol(inW)';
end

end
