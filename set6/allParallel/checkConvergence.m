close all
clear all


%% Check convergence
figure
Tn = sortrows(readtable('varyn.txt'))
ns = Tn{:,1};
loglog(ns, Tn{:,4})
hold on
loglog(ns, 1./ns.^2)
loglog(ns, 1./ns)
legend('Our data','quadratic convergence','linear convergence') 
xlabel('problem size N')
ylabel('max error')
saveas(gcf,'../report/figures/checkConv.png')


%% Check how runtime scales on a single node with only MPI or OMP
figure
Tomp = sortrows(readtable('t1-12.txt'))
Tmpi = sortrows(readtable('p1-12.txt'))
plot(Tmpi{:,2},Tmpi{:,5}, '*-')
hold on
plot(Tomp{:,3},Tomp{:,5}, '*-')
n = Tmpi{1,1};
xlabel('Number of threads/processors')
ylabel('Runtime in seconds')
string = sprintf('Runtime for N = %d \n with t = 1 (MPI implementation)', n);
string2 = sprintf('Runtime for N = %d \n with p = 1 (OMP implementation)',n);
legend(string, string2)
saveas(gcf,'../report/figures/runtime_either_MPI_OMP.png')

%% Check how many nodes is best
Tpn = sortrows(readtable('varynp.txt'),6)
figure
plot(Tpn{:,6},Tpn{:,5}, '*-')
xlabel('Number of nodes')
ylabel('Runtime in seconds')
string = sprintf('Runtime for N = %d \n with p = 12 and t = 1', Tpn{1,1});
legend(string)
saveas(gcf,'../report/figures/bestnodes.png')


%% p*t = 36
figure
Tt = sortrows(readtable('varyt.txt'))
plot(Tt{:,2},Tt{:,5}, '*-')
n = Tt{1,1};
xlabel('Number of MPI processors')
ylabel('Runtime in seconds')
string = sprintf('Runtime for N = %d \n with p*t = 36', n);
legend(string)
saveas(gcf,'../report/figures/pt36.png')




%% Task d - Report speedup
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



%% Check how runtime scales with n
figure
Tn = sortrows(readtable('varyn.txt'))
ns = Tn{4:end,1};
loglog(ns, Tn{4:end,5})
hold on
loglog(ns, ns.^2.*log(ns)/(ns(1)^2*log(ns(1))))
%loglog(ns, ns.^2/ns(1)^2)
%loglog(ns, ns/ns(1))
ylabel('Runtime in seconds')
xlabel('Problem size N')
string2 = sprintf('Runtime with p = %d and t = %d', Tn{1,2}, Tn{1,3});
legend(string2, 'N^2*log(N)')%, 'quadratic scaling','linear scaling');
saveas(gcf,'../report/figures/runtimeN.png')
