function [y] = testF3(x)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
y = 0;
for w = 1:1:20
    y = y + (1 - (w / 30) * (w / 30)) * sin(w * x) / w;
end
end