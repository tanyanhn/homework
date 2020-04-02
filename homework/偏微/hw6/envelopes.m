clear

syms x;
y1(x) = 1/x;
y2(x) = - 1/x^2;
%y3(x,i,k) = k(x - i) + y1(i);
fplot(y1,[-1,1]);
hold on;
for i = -1:0.11:1
k = y2(i);
%     y3(x,i) = k(x - i) + y1(i);
x = -1:0.11:1;
y = k.*(x - i) + y1(i);
plot(x,y);
hold on;
end
axis([-1 1 -1 1]);
