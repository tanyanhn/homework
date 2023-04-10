function [y] = testF2(x)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
y = 0;
for w = 1:1:20
    y = y + sin(w * x) / w;
end
end