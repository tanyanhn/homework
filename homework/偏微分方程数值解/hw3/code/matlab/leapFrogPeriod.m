function [vn] = leapFrogPeriod(v0, h, k, n, a)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
row = length(v0);
E = eye(row);

Q = ((k) / (h)) * circshift(E, [-1, 0]) - ((k) / (h)) * circshift(E, [1, 0]);

vn(:, 1:2) = v0; 
for i = 2:1:n
    vn(:, i + 1) = (Q * vn(:, i) + E * (1 - k * a) * vn(:, i - 1)) / (1 + k * a);  
end
end