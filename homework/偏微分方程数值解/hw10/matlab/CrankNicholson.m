function [vn] = CrankNicholson(v0, a, D0, DpDm, k, n)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
len = length(D0);
E = eye(len);
% E u - k / 2 * (D0 a * D0 u + E a * D+D_ u)
Al = E - k / 2 * D0 * a .* D0 - k / 2 * E * a .* DpDm;
Ar = E + k / 2 * D0 * a .* D0 + k / 2 * E * a .* DpDm;
A = Ar / Al;
vn(:, 1) = v0; 
for i = 1:1:n
    vn(:, i + 1) = Al \ (Ar * vn(:, i));  
end
end