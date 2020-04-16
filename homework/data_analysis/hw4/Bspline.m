clear;
clc;

syms x;
y0 = 0;
y1 = (x+1)*(x+1)/2;
y2 = (1+2*x-2*x^2)/2;
y3 = (2-x)^2/2;

fplot(y0,[-2,-1]);
hold on;
fplot(y1,[-1,0]);
hold on;
fplot(y2,[0,1]);
hold on;
fplot(y3,[1,2]);
hold on;
fplot(y0,[2,3]);
