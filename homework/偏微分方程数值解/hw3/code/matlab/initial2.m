function [v0] = initial2(f, h, k, domain, a)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
length = domain(2) - domain(1);
row = single(length / h);
v0 = zeros(row + 1, 2);

for i = 0:1:row
    v0(i + 1, 1) = f(i * h + domain(1));
end

E = eye(row + 1);

Q = E + (k / (2 * h)) * circshift(E, [-1, 0]) - ( k / (2 * h)) * circshift(E, [1, 0]);
v0(:, 2) = Q * v0(:, 1) - k * a .* v0(:, 1);
end