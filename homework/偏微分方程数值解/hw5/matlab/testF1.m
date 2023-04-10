function [y] = testF1(x)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
x  = mod(x, 2 * pi);
y = (pi - x) / 2;
end