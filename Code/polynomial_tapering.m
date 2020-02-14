function [y, dydx] = polynomial_tapering(x, x_lower,x_upper)
% function [y, dydx] = polynomial_tapering(x, x_lower,x_upper)
% Let f(x) = a*x^3 + b*x^2 + c*x + d, and satisfy
% f(x_lower) = 0; f(x_upper) = 1; f'(x_lower) = f'(x_upper) = 0

a = 2 ./ (x_lower.^3 - x_upper.^3 + 3 * (x_lower) .* (x_upper.^2) - 3 * (x_lower.^2) .* (x_upper));
b = -1.5 * a .* (x_lower + x_upper);
c = 3 * a .* x_lower .* x_upper;
d = 0.5 * a .* (x_lower.^3) - 1.5 * a .* (x_lower.^2) .* (x_upper);


y = a .* (x.^3) + b .* (x.^2) + c .* x + d;

dydx = 3 * a .* (x.^2) + 2 * b .* x + c;



