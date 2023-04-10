function [output] = Dzero(A)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

output = (Dplus(A) + Dminus(A)) / 2;

end