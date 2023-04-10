function [ret] = L2(A)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
ret = sum(A .* A);
row = size(A);
ret = sqrt(ret / row(1));
end