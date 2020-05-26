clc;
clear;
clear all;

syms x;

y1 = 1 / x * exp(-x);
y2 = x * exp(-x) / (1 - exp(-x));

fplot(x,y1,[0.1,1]);
hold on;
fplot(x,y2,[0,1]);
