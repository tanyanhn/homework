function [vn] = laxWendroffPeriod(v0, h, k, n)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
row = length(v0);
E = eye(row);

Q = (1 - (k * k) / (h * h)) * E ...
    + ((k * k) / (2 * h * h) + k / (2 * h)) * circshift(E, [1, 0])...
    + ((k * k) / (2 * h * h) - k / (2 * h)) * circshift(E, [-1, 0]);

vn =  Q^n * v0;  
end