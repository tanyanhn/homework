% Program to plot stability plots for Runge Kutta methods
% For details on plotting stability regions for One-step Method,
% see pg 163 of Randall J. LeVeque book on " Finite Difference
% Methods for Ordinary and Partial Differential Equations


% In case on one-step method, any time-step n+1 can be written as
% u(n+1) = R(z)*u(n)
% where R(z) is some function of z = lambda*dt

% The region of absolute stability for a one-step method is simply the
% region where |R(z)| <= 1


% We demonstrate here stability curve for RK2, RK3 and RK4 methods.

% Setting up the grid points
clear
axisbox = [-4 4 -4 4];
xa = axisbox(1); xb = axisbox(2); ya = axisbox(3); yb = axisbox(4);
nptsx = 51;
nptsy = 51;
x = linspace(xa,xb,nptsx);
y = linspace(ya,yb,nptsy);
[X,Y] = meshgrid(x,y);
Z = X + 1i*Y;

% 2nd order Runge-Kutta method
RK2 = @(z) (1 + z + 1/2*z.^2);

% 3rd order Runge-Kutta method
RK3 = @(z) (1 + z + 1/2*z.^2 + 1/6*z.^3);

% 4th order Runge-Kutta method
RK4 = @(z) (1 + z + 1/2*z.^2 + 1/6*z.^3 + 1/24*z.^4);

% Evaluation of R(z) inside axisbox
Rval2 = RK2(Z); Rval3 = RK3(Z); Rval4 = RK4(Z);
%----------------------------------------------------------------------
% evaluation of |R(z)| for finding S, region of absolute stability:

Rabs2 = abs(Rval2); Rabs3 = abs(Rval3); Rabs4 = abs(Rval4);

% Plotting the stability curve using the contour plot of MATLAB

figure()
clf()
hold on
% plot axes:
plot([xa xb],[0 0],'k-','LineWidth',2)
plot([0 0],[ya yb],'k-','LineWidth',2)
contour(x,y,Rabs2,[1 1],'r-','LineWidth',2);
contour(x,y,Rabs3,[1 1],'b-','LineWidth',2);
contour(x,y,Rabs4,[1 1],'g-','LineWidth',2);
legend('RK2','RK3','RK4')
title('Region of absolute stability','FontSize',15)
axis([xa xb ya yb])
set(gca,'FontSize',15)
hold off

% store data
c2 = contourc(x,y,Rabs2,[1 1])';
c2(1,:) = [];
c3 = contourc(x,y,Rabs3,[1 1])';
c3(1,:) = [];
c4 = contourc(x,y,Rabs4,[1 1])';
c4(1,:) = [];
save "RK_stability.dat" c2 c3 c4