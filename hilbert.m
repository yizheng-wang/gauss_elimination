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
xlabel('Hilbert����')
ylabel('�������µľ������')
title('Hilbert�������������Ԫ���͸�˹��ȥ���ľ������Ĺ�ϵ')