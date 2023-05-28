function [output] = Dminus(A)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

output = A - circshift(A, [1, 0]);

end