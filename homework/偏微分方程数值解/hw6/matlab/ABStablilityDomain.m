
clear


f = @(z, u)(1 - 6/11 * u) * z^3 - 18 /11 * z^2 + 9 /11 * z - 2 /11;

dist = 0.1;
bound = 4;
lowbound = 0.01;
axisbox = [-bound bound -bound bound];
xa = axisbox(1); xb = axisbox(2); ya = axisbox(3); yb = axisbox(4);
nptsx = bound / dist * 2 + 1;
nptsy = nptsx;
xl = linspace(xa,xb,nptsx);
yl = linspace(ya,yb,nptsy);
[X,Y] = meshgrid(xl,yl);
uf = zeros(size(X));

i = 0;
j = 0;
for xu = -bound:dist:bound
    i = i + 1;
    j = 0;
    for yu = -bound:dist:bound
        j = j + 1;
        zu = xu + 1i * yu;
        for x = -bound:dist:bound
            for y = -bound:dist:bound
                z= x + 1i * y;
                if (abs(f(z, zu)) < lowbound)
                    uf(i, j) = max(uf(i, j), abs(zu));
                end
    
            end
        end

    end
end

plot([xa xb],[0 0],'k-','LineWidth',2)
plot([0 0],[ya yb],'k-','LineWidth',2)
contour(xl,yl,uf,[1 1],'r-','LineWidth',2);
