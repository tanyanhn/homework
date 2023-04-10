function FTCS(h_0,k_0,Total_Time_0,a,eta,Initial_Condition,str_name_type)
% FUNCTION - FTCS scheme for $u_t + a * u_x = \eta * u_{xx}$.
% INPUTS:
%   h_0 - Double: Space step
%   k_0 - Double: Time step
%   Total_Time_0 - Double: The end time for time advance
%   a - Double: Used in $u_t + a * u_x = \eta * u_{xx}$
%   Intial_Condition - Function Handle: The initial condition
%   str_name_type - String: description of the plot

Time_steps = round(Total_Time_0/k_0);
N = round(2*pi/h_0)-1;
h = h_0;
k = k_0;
u = zeros(1,N+1);
alpha = 2*eta*k/(h*h);
lambda = a * k / h;

% Set Initial Condition
for i = 1 : (N+1)
    u(i) = Initial_Condition(i*h);
end

% Time Advance by FTCS
for i= 1 : Time_steps
    u_temp=zeros(1,N+1);
    for j = 2 : N
        u_temp(j) = u(j) + alpha/2*(u(j+1)+u(j-1)-2*u(j)) - lambda/2*(u(j+1)-u(j-1));
    end
    u_temp(1) = u(1) + alpha/2*(u(2)+u(N+1)-2*u(1)) - lambda/2*(u(2)-u(N+1));
    u_temp(N+1) = u(N+1) + alpha/2*(u(1)+u(N)-2*u(N+1)) - lambda/2*(u(1)-u(N));
    u=u_temp;
end

str_title = sprintf('(N,\\alpha,\\lambda)= (%d,%.2f,%.2f)',N,alpha,lambda);
figure(1);
axis([h 2*pi 1.05*min(u) 1.05*max(u)]);
hold on;
Space_Grid = (1 : N+1) *h;
fp1=plot(Space_Grid,u,'k');
fp1.LineStyle='-';
fp1.LineWidth = 0.4;
hold on;
str_postfix_png='.png';
title(str_title);
saveas(1,[sprintf('eta=%.2f:N=%d,',eta,N),str_name_type,str_postfix_png]);
clf(1);
end