clc;
clear;

file = "output/output";

syms x;
g = 1 / (1 + 25 * x^2);


for i=1:1:5
    figure;
    s = "N=";
    s = s + i;
    title("s");
    filename = file + i;
    data = dlmread(filename);
    [m,n] = size(data);
%    error(i,1) = 0;
    for j=1:1:m/3
       a = data(3 * j - 2,1);
       b = data(3 * j - 2,2);
       x0 = data(3 * j - 1,1);
       param = data(3 * j,:);
      % syms x;
       y = param(1,1) * (x - x0)^3 + param(1,2) * (x - x0)^2 + param(1,3) * (x - x0)^1 + param(1,4);
       fplot(y,[a,b]);
       hold on;
       e = g - y;
       fplot(e,[a,b]);
%        k =  max(subs(e,x,(a:0.001:b)));
%        if(error(i,1) < k)
%            error(i,1) = k;
%        end
       hold on;
    end
    %gtext('interpolation curve');
    %gtext('error curve');
end

