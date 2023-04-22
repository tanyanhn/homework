
clear 
close all

pi2 = 2 * pi;
domain = [0, pi2];
time = pi;
N = 32;
% k = 0.001;
a = 1;
eta = 0.1;

alpha = 1;
lambda = 1;
ample = 32;
exactSolution = exactSol(@testF1, time, N * ample, domain, a, eta);

for ifunc = 1:3

N = N * 2;
ample  = ample / 2;
h = pi2 / N;

x = domain(1):h:domain(2);
k1 = alpha * h * h / (2 * eta);
k2 = lambda * h / a;
k = min(k1, k2);
k = min(k , 0.001);
n = single(time / k);

    for p = 1:2
        Qc = convection(h, N, 1);
        Qd = diffusion(h, N, p);
        Q = -a * Qc + eta * Qd;
        v0 = initial2(@testF1, Q, h, k, domain);
        vn = oneStepPeriod(v0, k, Q, n);
        figure

        % fig = plot(x, (vn(:, end)),'-'); hold on
        fig = plot(x, (vn(:, end) - exactSolution(1:ample:end)),'-'); hold on
        strtitle  = sprintf('grid-%d-Q_%d-result',N , 2 * p);
        title(strtitle);
        saveas(fig, ['../fig/',strtitle, '.png']);
    end
end

