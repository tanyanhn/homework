function [v0] = initial2(f, Q, h, k, domain)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
length = domain(2) - domain(1);
row = single(length / h);
v0 = zeros(row + 1, 2);

for i = 0:1:row
    v0(i + 1, 1) = f(i * h + domain(1));
end

v0(:, 2) = k * Q * v0(:, 1) + v0(:, 1);
end