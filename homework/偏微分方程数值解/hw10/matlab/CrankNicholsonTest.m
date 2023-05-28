
clear 
close all

pi2 = 2 * pi;
domain = [0, pi2];
time = pi / 100;
N = 512;
h = pi2 / N;
k = h / 10;
x = 0:h:pi2;
E = eye(N + 1);

% u_t = a_x u_x + a u_{xx} 
% u_t = D0 a * D0 u + a * D+D_ u

D0 = Dzero(E) * h^-1;
DpDm = Dplus(Dminus(E)) * h^-2;

n = single(time / k);

testF = {@testF1};
testa = {@ax, @ax2};
extactSol ={@sol1, @sol2};
for ifunc = 1:1
    for iax = 1:2
        v0 = initial1(testF{ifunc}, h, domain);
        a = initial1(testa{iax}, h, domain);
        vn = CrankNicholson(v0, a, D0, DpDm, k, n);
        figure
        fig = plot(x, vn(:, end),'-'); hold on
        strtitle  = sprintf('a(x)_%d_calculateResult',iax);
        title(strtitle);
        saveas(fig, ['../fig/',strtitle, '.png']);
        figure 
        fig = plot(x, extactSol{iax}(x, time));
        strtitle  = sprintf('a(x)_%d_exactSolution',iax);
        title(strtitle);
        saveas(fig, ['../fig/',strtitle, '.png']);
        max(vn(:, end))
    end
end

