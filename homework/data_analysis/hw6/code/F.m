clc;
clear;

file = "output/F";

for i=1:1:2
    %figure;
    s = "N=";
    s = s + i;
    title("s");
    filename = file + i;
    data = dlmread(filename);
    [m,n] = size(data);
    subplot(1,2,i);
    plot(data(1,:),0,'+b');
    if(i ==2)
        hold on;
        plot(data(3,:),0,'*r');
    end
    %hold on;
end