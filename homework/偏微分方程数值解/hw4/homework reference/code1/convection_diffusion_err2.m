a = -1;
eta = 0.5;
x_max = 2*pi;
N = 10;
h = x_max/N;
k = 0.5*h*h/(2*eta);
mu = a*k/h;
max_iter = round(-N/mu);
numfine = 3;

hold on;
for i = 1:numfine+1
    u = initial_data(N,x_max);
    res = convection_diffusion(u,N,mu,x_max,a,eta,max_iter);
    plot_fig(x_max,N,res);
    N = N*2;
    h = x_max/N;
    k = 0.5*h*h/(2*eta);
    mu = a*k/h;
    max_iter = round(-N/mu);
end
legend('N = 10','N = 20','N = 40','N = 80');