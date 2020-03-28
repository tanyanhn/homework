clc
clear

syms x;
y = 1 /(1 + x*x);
fplot(y,[-5,5])
gtext('orignal line');
hold on;

for i = 2:2:8
    s = "Data/newton " + i;
    d = load(s);
    parameter1 = d(1,:);
    parameter2 = d(2,:);
    x = -5:0.001:5;
    y = f(parameter1,parameter2,x);    
    plot(x,y) 
    text = [num2str(i),'th'];
    gtext(text);
    hold on;
end

title('Newton interpolation');
hold off;

clear

figure;
syms x;
y = 1 /(1 + 25*x*x);
fplot(y,[-1,1])
gtext('orignal line');
hold on;

for i = 5:5:20
    s2 = "Data/newton2 " + i;
    d2 = load(s2);
    parameter1 = d2(1,:);
    parameter2 = d2(2,:);
    x = -1:0.001:1;
    y2 = f(parameter1,parameter2,x);    
    plot(x,y2) 
    text = [num2str(i),'th'];
    gtext(text);
    hold on;
end

title('Newton interpolation 2');
hold off;

figure;
syms x;
y = 1 /(1 + 25*x*x);
fplot(y,[-1,1])
gtext('orignal line');
hold on;

for i = 5:5:20
    s2 = "Data/chebyshev3 " + i;
    d2 = load(s2);
    parameter1 = d2(1,:);
    parameter2 = d2(2,:);
    x = -1:0.001:1;
    y2 = f(parameter1,parameter2,x);    
    plot(x,y2) 
    text = [num2str(i),'th'];
    gtext(text);
    hold on;
end

title('Chebyshev interpolation 3');
hold off;

