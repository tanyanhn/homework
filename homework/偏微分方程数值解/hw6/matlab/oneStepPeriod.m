function [vn] = oneStepPeriod(v0, k, Q, n)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% len = length(Q);
% E = eye(len);
Q = k * Q;
vn(:, 1:2) = v0; 
for i = 1:1:n
    vn(:, i + 1) =  Q * vn(:, i) + vn(:, i);  
end
end