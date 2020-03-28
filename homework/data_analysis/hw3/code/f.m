function [a] = f(parameter1,parameter2,x)
a = 0;
i = length(parameter1(:));
while i > 0
    a = parameter1(i) + a .* (x - parameter2(i));
    i=i-1;
end


