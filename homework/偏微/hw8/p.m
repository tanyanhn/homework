clear;
clc;
close;

figure;

syms x;
syms y;
y1 = 2*x + 2;
fplot(y1,[-1,0]);
hold on;
y2 = 1/2*x;
fplot(y2,[0,2]);
hold on;
y3 = x - 1;
fplot(y3,[1,2]);
hold on;
y4 = 1/4 * x^2;
fplot(y4,[2, 4 + 2* sqrt(2)]);
hold on;
y5 = 2*x - 2;
fplot(y5,[4 + 2 * sqrt(2), 9]);
hold on;
z = 2:0.01:(6 + 4 * sqrt(2));
y6 = z;
x1 = z - sqrt(2*z);
plot(x1,y6);
hold on;

arrow([-2,0],[0,2]);
arrow([-2 - 2* sqrt(2),0],[4 + 2*sqrt(2),6 + 4*sqrt(2)]);
arrow([0,0],[4 + 2*sqrt(2),6 + 4*sqrt(2)]);
arrow([4 + 2*sqrt(2),0],[4 + 2*sqrt(2),6 + 4*sqrt(2)]);
arrow([-6,0],[8,14]);
arrow([8,0],[8,14]);
a = 4;
b = a - sqrt(2*a);
arrow([b-4,0],[b,a]);
arrow([0,0],[b,a]);
a = 7;
b = a - sqrt(2*a);
arrow([b-7,0],[b,a]);
arrow([0,0],[b,a]);
arrow([-0.5,0],[-0.5,1]);
arrow([0,0],[0,2]);
arrow([0.5,0],[3/2,0.5]);
arrow([3/2,0],[3/2,0.5]);
a= 3;
b = 1/4* a*a;
arrow([a,0],[a,b]);
arrow([0,0],[a,b]);
a= 4.2;
b = 1/4* a*a;
arrow([a,0],[a,b]);
arrow([0,0],[a,b]);
a= 5.5;
b = 1/4* a*a;
arrow([a,0],[a,b]);
arrow([0,0],[a,b]);







