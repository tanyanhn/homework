function u = fill_ghostcell_period(u,N,num)
count = 0;
while count < num
    u = [u(N),u,u(2+2*count)];
    count = count + 1;
end
    
    
