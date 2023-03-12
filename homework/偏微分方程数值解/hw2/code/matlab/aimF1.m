function [y] = aimF1(x, t)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
x  = mod(x, 2 * pi);
y = x.*(x < pi) + (2 * pi - x) .* (x > pi);
end