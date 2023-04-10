function [y] = extSolution(f, x, t, a)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
x = mod(x, 2 * pi);
y = exp(-a * t) .* f(x + t);
end