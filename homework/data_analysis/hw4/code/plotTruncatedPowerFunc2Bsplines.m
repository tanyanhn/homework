function plotTruncatedPowerFunc2Bsplines(n)
   % f = zero(n+2,n+2);
   x = -1:0.01:(n+3);
   y = zeros(n+2,n+2,(n+4)*100 +1);
    for j = 1:1:(n+2)
        for i = j:1:(n+2)
            if j ==1
                y(i,j,:) = truncate(x,n,i);
                for k = 1:1:((n+4)*100 + 1)
                    m(k) = y(i,j,k);
                end
                %m = y(i,j,:);
                figure;
                plot(x,m);
                axis([0 n+3 0 1]);
                %plot(x, y(i,j,:));
            else
                y(i,j,:) = (y(i,j-1,:) - y(i-1,j-1,:))./(j-1) ;
                for k = 1:1:((n+4)*100 + 1)
                    m(k) = y(i,j,k);
                end
                figure;
                plot(x,m);
                axis([0 n+3 0 1]);
            end
        end
    end
    figure;
    %plot(n+1.*m);
    %axis([0 n+3 0 1]);
end

function y = truncate(x,k,x0)
    %x = (x0-k-2):0.01:(x0 +1);
    y = ((x0 - x).^k).*(x < x0) + (0).*(x >= x0);
    %fun = @c;
end