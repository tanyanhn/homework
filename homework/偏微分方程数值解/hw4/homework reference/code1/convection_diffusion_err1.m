x_max = 2*pi;
N = 10;
mu = -0.5;
max_iter = round(-N/mu);
numfine = 1;
a = -1;
eta = 0.5;

hold on;
for i = 1:numfine+1
    u = initial_data(N,x_max);
    res = convection_diffusion(u,N,mu,x_max,a,eta,max_iter);
    plot_fig(x_max,N,res);
    N = N*2;
    max_iter = max_iter*2;
end
legend('N = 10','N = 20');
u = initial_data(N,x_max);
res = convection_diffusion(u,N,mu,x_max,a,eta,max_iter)