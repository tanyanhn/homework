function [y] = testF2(x)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
x  = mod(x, 2 * pi);
y = sin(x) +sin(10 * x);
end