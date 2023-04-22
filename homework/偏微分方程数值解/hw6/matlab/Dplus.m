function [output] = Dplus(A)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

output = circshift(A, [-1, 0]) - A;

end