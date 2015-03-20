close all
clear all

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
Tt = sortrows(readtable('varyt.txt'))
plot(Tt{:,2},Tt{:,5})
n = Tt{1,1};
xlabel('Number of processors')
ylabel('Runtime in seconds')
string = sprintf('Runtime for n = %d \n with p*t = 36', n);
legend(string)

figure
plot(Tn{:,1}, Tn{:,5})
ylabel('Runtime in seconds')
xlabel('Problem size n')
string2 = sprintf('Runtime with p = %d and t = %d', Tn{1,2}, Tn{1,3});
legend(string2);

% Task d - Report speedup