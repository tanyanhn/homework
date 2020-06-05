clc;
clear;

for i = 1:1:2
    %figure;
    file = "data/matrix";
    %s = s + i;
    %title("s");
    filename = file + i;
    data = dlmread(filename);
    %[m,n] = size(data)
    a = zeros(3,3);
    for j = 1:1:3
        for k = 1:1:3
            a(j,k) = data(1,j + (k-1) * 3);
        end
    end
    a
    b = inv(a)
    x = norm(a,2)
    y = norm(b,2)
    x * y
end