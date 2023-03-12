
% clear 
% close all
% 
% pi2 = 2 * pi;
% domain = [0, pi2];
% time = pi2;
% 
% exactSolution;

% forwardEulerPeriod with h = 2*pi/10
clear h k n v0 t
h = pi2 / 10;
k = h / 2;
n = single(time / k);
v0 = initial(@testF1, h, domain);
vn10 = forwardEulerPeriod(v0, h, k, n);
j = jump;
lm = max(abs((vn10 - realvn(1:j:end, 1:j:end))));
figure
t = (0:1:length(lm) - 1) * k;
plot(t, lm,'-');

% h = 2*pi/100
clear h k n v0 t
h = pi2 / 100;
k = h / 2;
n = single(time / k);
v0 = initial(@testF1, h, domain);
vn100 = forwardEulerPeriod(v0, h, k, n);
j = jump / 10;
lm = max(abs((vn100(1:1:end, 1:1:end) - realvn(1:j:end, 1:j:end))));
figure
t = (0:1:length(lm) - 1) * k;
plot(t, lm,'-');