close all
clear all

errors = [1.927711e-01 , 4.476185e-02, 1.098931e-02 , 2.734955e-03 , 6.829684e-04 , 1.706940e-04 , 4.267049e-05, 1.066744e-05, 2.666847e-06 ];
ns = [4, 8, 16, 32, 64, 128, 256, 512, 1024];

Tn = sortrows(readtable('varyn.txt'))
ns = Tn{:,1};
errors = Tn{:,4};

loglog(ns, errors)
hold on
loglog(ns, 1./ns.^2)
loglog(ns, 1./ns)
legend('Our data','quadratic convergence','linear convergence') 
xlabel('problem size n')
ylabel('max error')


figure
Tp = sortrows(readtable('varyp.txt'))
plot(Tp{:,2},Tp{:,5})
n = Tp{1,1};
xlabel('Number of processors')
ylabel('Runtime in seconds')
string = sprintf('Runtime for n = %d \n with p*t = 36', n);
legend(string)

