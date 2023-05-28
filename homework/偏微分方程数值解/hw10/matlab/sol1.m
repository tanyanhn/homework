function [y] = sol1(x, t)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
x  = mod(x, 2 * pi);
y = exp(-t) * sin(x);
end