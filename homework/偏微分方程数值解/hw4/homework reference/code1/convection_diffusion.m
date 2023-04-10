function u = convection_diffusion(u,N,mu,x_max,a,eta,num)
for i = 1:num
    u = convection_diffusion_onestep(u,N,mu,x_max,a,eta);
end