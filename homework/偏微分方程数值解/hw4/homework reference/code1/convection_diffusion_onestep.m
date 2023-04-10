function res = convection_diffusion_onestep(u,N,mu,x_max,a,eta)
h = x_max/N;
k = mu*h/a;
u = fill_ghostcell_period(u,N,1);
tmp = k*eta/(h*h);
res = (tmp+0.5*mu)*u(1:N+1) + (1-2*tmp)*u(2:N+2) + (tmp-0.5*mu)*u(3:N+3);