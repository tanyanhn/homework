function [vn, alpha, lambda] = scheme(v0, h, k, timeStep, a, eta)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
row = length(v0);
E = eye(row);

alpha2 = (k * eta) / (h * h);
lambda2 = (a * k) / (2 * h);

alpha = alpha2 * 2;
lambda = lambda2 * 2;

Q = E + alpha2 * (circshift(E, [1, 0]) - 2 * E + circshift(E, [-1, 0])) ...
    - lambda2 * (circshift(E, [1, 0]) - circshift(E, [-1, 0]));

vn(:, 1) = v0;
for i = 2:(timeStep + 1)
    vn(:, i) = Q * vn(:, i - 1);  
end
end