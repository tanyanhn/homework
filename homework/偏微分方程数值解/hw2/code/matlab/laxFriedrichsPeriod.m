function [vn] =  laxFriedrichsPeriod(v0 ,h, k, n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
row = length(v0);
E = eye(row);

Q = ((h + k) / (2 * h)) * circshift(E, [1, 0]) + ((h - k) / (2 * h)) * circshift(E, [-1, 0]) ;

vn(:, 1) = v0;
for i = 2:(n + 1)
    vn(:, i) = Q * vn(:, i - 1);  
end
end