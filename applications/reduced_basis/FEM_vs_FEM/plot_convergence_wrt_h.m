clear all
clc

load('errors_convergence_wrt_h.mat')
figure(5)
semilogy(1:length(errs),errs,'.-','Markersize',10)

legend('h = 1/20','h = 1/40','h = 1/80','h = 1/160');

xlabel('number of modes')
ylabel('L2 error')

grid on
axis square
axis([1 16 1e-10 1e-3])
set(gca,'fontsize',15)

title('Conformal mesh')