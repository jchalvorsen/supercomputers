close all
clear all

Tn = sortrows(readtable('varyn.txt'))
ns = Tn{:,1};
loglog(ns, Tn{:,4})
hold on
loglog(ns, 1./ns.^2)
loglog(ns, 1./ns)
legend('Our data','quadratic convergence','linear convergence') 
xlabel('problem size n')
ylabel('max error')
saveas(gcf,'../report/figures/checkConv.png')

% Check how runtime scales with n
figure
loglog(ns, Tn{:,5})
hold on
loglog(ns, ns.^2)
loglog(ns, ns)
ylabel('Runtime in seconds')
xlabel('Problem size n')
string2 = sprintf('Runtime with p = %d and t = %d', Tn{1,2}, Tn{1,3});
legend(string2, 'quadratic scaling','linear scaling');

figure
Tt = sortrows(readtable('varyt.txt'))
plot(Tt{:,2},Tt{:,5})
n = Tt{1,1};
xlabel('Number of processors')
ylabel('Runtime in seconds')
string = sprintf('Runtime for n = %d \n with p*t = 36', n);
saveas(gcf,'../report/figures/runtime.png')
legend(string)



% Task d - Report speedup
figure
Tp = sortrows(readtable('varyp.txt'))
p = Tp{:,2};
t = Tp{:,5};
Sp = t(1)./t;
plot(p, Sp, '*-')
ylabel('Speedup')
xlabel('Processors')
saveas(gcf,'../report/figures/speedup.png')

figure
np = Sp./p;
plot(p, np, '*-')
ylabel('Parallel efficiacy');
xlabel('Processors');
saveas(gcf,'../report/figures/efficiacy.png')
