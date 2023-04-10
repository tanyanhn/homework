function [y] = testF3(x)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
x  = mod(x, 2 * pi);
y = sin(x) +sin(100 * x);
end