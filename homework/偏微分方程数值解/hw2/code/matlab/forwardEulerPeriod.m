function [vn] = forwardEulerPeriod(v0, h, k, n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
row = length(v0);
E = eye(row);

Q = E + (k / (2 * h)) * circshift(E, [1, 0]) + (- k / (2 * h)) * circshift(E, [-1, 0]) ;
vn(:, 1) = v0;
vi = 2:1:(n + 1);
for i = vi
    vn(:, i) = Q * vn(:, i - 1);  
end

size(vn)