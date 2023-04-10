
clear 
close all

pi2 = 2 * pi;
domain = [0, pi2];
time = pi;
N = 128;
k = 0.001;
h = pi2 / N;
x = 0:h:pi2;
E = eye(N + 1);



n = single(time / k);
Q{1}(:, :) = Dzero(E) * h^-1;
Q{2}(:, :) = Q{1}(:, :) - h^-1 / 6 * Dzero(Dplus(Dminus(E)));
Q{3}(:, :) = Q{2}(:, :) + h^-1 / 30 .* Dzero(Dplus(Dplus(Dminus(Dminus(E)))));

testF = {@testF1, @testF2, @testF3};
for ifunc = 1:3
    for p = 1:3
        v0 = initial2(testF{ifunc}, -Q{p}, h, k, domain);
        vn = leapFrogPeriod(v0, -Q{p}, h, k, n);
        figure
        fig = plot(x, vn(:, end),'-'); hold on
        strtitle  = sprintf('initial_function_%d-Q%d_result',ifunc, 2 * p);
        title(strtitle);
        saveas(fig, ['../fig/',strtitle, '.png']);
    end
end

