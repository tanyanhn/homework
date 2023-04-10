% a = 1; TotalTime = pi;
% Initial Condition y = sin(x0)+sin(4*x0)
clear
Total_Time = pi;
a = 1;
alpha_max = 1;
lambda_max = 1;
for eta = [1.0,0.01]
    for N = [10,20,40]
        h = 2*pi/(N+1);
        %k = min(alpha_max*h*h/(2*eta),abs(lambda_max*h/a));
        k1 = alpha_max*h*h/(2*eta);
        k2 = lambda_max*h/a;
        FTCS(h,k1,Total_Time,a,eta,@LowFrequency,'alpha=1');
        FTCS(h,k2,Total_Time,a,eta,@LowFrequency,'lambda=1');
    end
end


function y = LowFrequency(x0)
y = sin(x0)+sin(4*x0);
end
