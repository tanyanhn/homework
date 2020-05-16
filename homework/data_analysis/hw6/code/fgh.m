clc;
clear;

file = "output/ma";

i = 0.02/100;
x = 0.99 : i : 1.01;

for i=1:1:3
    %figure;
    s = "N=";
    s = s + i;
    title("s");
    filename = file + i;
    data = dlmread(filename);
    [m,n] = size(data);
    subplot(1,3,i);
    plot(x, data);
    %hold on;
end