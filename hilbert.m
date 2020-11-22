for n = 1:40

H = hilb(n);
xex = ones(n, 1);
b = H * xex;
xauto = gauss(H,b);
xmax = gauss(H,b,1);
error_xauto(n) = norm(xex-xauto);
error_xmax(n) = norm(xex-xmax);
end
plot(1:40, error_xauto,'r')
hold on 
plot(1:40, error_xmax)
legend('auto','max')
xlabel('Hilbert阶数')
ylabel('二范数下的绝对误差')
title('Hilbert矩阵阶数与列主元法和高斯消去法的绝对误差的关系')