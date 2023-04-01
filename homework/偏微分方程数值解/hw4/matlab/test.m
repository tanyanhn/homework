
clear 
close all
% 

pi2 = 2 * pi;
domain = [0, pi2];
time = pi2;
a = 1;
eta = 1;
% exactSolution;

alpha = 1;
lambda = 1;
i = 1;
for N = [10 20 40]
    clear h k1 k2 timeStep1 timeStep2 u0 t
    h = 2 * pi / N;
    k1 = alpha * h * h / (2 * eta);
    k2 = lambda * h / a;
    timeStep1 = single(time / k1);
    timeStep2 = single(time / k2);

    u0 = initial(@initialValue, h, domain);
    [u1n, alpha1, lambda1] = scheme(u0, h, k1, timeStep1, a, eta);
    [u2n, alpha2, lambda2] = scheme(u0, h, k2, timeStep2, a, eta);

    figure
    Space_Grid = (1 : N+1) *h;
    fig = plot(Space_Grid, u1n(:, end), 'k');
    str_title = sprintf('(N,\\alpha,\\lambda)= (%d,%.2f,%.2f)',N,alpha1,lambda1);
    title(str_title);
    name = "fig/alpha_bound_" + num2str(i);
    saveas(fig, name, 'png');


    figure
    Space_Grid = (1 : N+1) *h;
    fig = plot(Space_Grid, u2n(:, end), 'k');
    str_title = sprintf('(N,\\alpha,\\lambda)= (%d,%.2f,%.2f)',N,alpha2,lambda2);
    title(str_title);
    name = "fig/lambda_bound_" + num2str(i);
    saveas(fig, name, 'png');

    i = i + 1;
end
