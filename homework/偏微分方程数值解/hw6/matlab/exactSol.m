function [sol] = exactSol(f, time, N, domain, a, eta);

h = (domain(2) - domain(1)) / N;
x = domain(1):h:domain(2);


alpha = 1;
lambda = 1;
k1 = alpha * h * h / (4 * eta);
k2 = lambda * h / a;
k = min(k1, k2);
n = single(time / k);


Qc = convection(h, N, 1);
Qd = diffusion(h, N, 1);
Q = -a * Qc + eta * Qd;
v0 = initial2(f, Q, h, k, domain);
vn = oneStepPeriod(v0, k, Q, n);
figure
sol = vn(:, end);
fig = plot(x, vn(:, end),'-'); hold on
strtitle  = sprintf('exactSolution');
title(strtitle);
saveas(fig, ['../fig/',strtitle, '.png']);

end