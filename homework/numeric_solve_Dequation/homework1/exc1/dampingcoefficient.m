n = 64;
k = 1:(n-1);
w = 2/3;
lambda = 1 - 2 * w * sin(k * pi/ 2 / n).^2;
n2 = 1:(n/2);
ck = cos(n2 * pi / 2 / n).^2;
sk = sin(n2 * pi / 2 / n).^2;

v = [0 0; 0 2; 1 1; 2 0; 2 2; 4 0];
for iv = 1:6
    x1 = 1:(n/2);
    x2 = (n/2):(n-1);
    c1 = lambda(x1) .^(v(iv,1) + v(iv,2)) .* sk;
    c2 = lambda(n - x2) .^v(iv,1) .* lambda(x2) .^v(iv,2) .* sk(n - x2);
    c3 = lambda(n - x1) .^v(iv,1) .* lambda(x1) .^v(iv,2) .* ck;
    c4 = lambda(x2) .^(v(iv,1) + v(iv,2)) .* ck(n - x2);
    subplot(3,2,iv)
    ax1=plot(x1,c1,'r'); hold on;
    ax2=plot(x2,c2,'g'); hold on;
    ax3=plot(x1,c3,'b'); hold on;
    ax4=plot(x2,c4,'m'); hold on;
    legend([ax1,ax2,ax3,ax4],{'c1','c2','c3','c4'});
    title("v1 = " + num2str(v(iv,1)) + ", v2 = " + num2str(v(iv,2))); 
    xlabel('k');
end
% w = zeros(n-1,n+1);
% for k = 1:(n-1)
%     lambda(k) = cos(k * pi/n) * cos(k * pi/n);
%     for j = 0:n
%         w(k,j+1) = sin(j * k * pi / n);
%         for tj = 1:j
%             w(k,j+1) = w(k,j+1) * cos(k * pi /n);
%         end
%     end
% end

