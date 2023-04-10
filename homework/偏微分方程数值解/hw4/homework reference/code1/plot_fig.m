function plot_fig(x_max,N,res)
h = x_max/N;
m = size(res,2);
plot(0:h:(m-1)*h,res,'.-');xlim([0 2*pi]);ylim([-0.2 0.2]);