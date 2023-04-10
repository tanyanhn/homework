function u = initial_data(N,x_max)
u = [];
for x = 0:x_max/N:x_max
    u = [u,initial_function(x)];
end