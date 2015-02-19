#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

double * generateV(n){
    // Generate the array of data:
    double *v = malloc(sizeof(double) * n);
    int i;
    for (i = 1; i <= n; ++i){
        v[i-1] = 1/(i*i*1.0);
    }
    return v;
}

double computeSum(double *v, int n){
    double sum=0.0;
    int i;

    #pragma omp parallel for schedule(static) reduction(+:sum)
    for (i = 0; i < n; ++i){
        sum += v[i];
    }
    return sum;
}

int main(int argc, char **argv)
{
    int k, n;
    double S = M_PI*M_PI/6;
    double sum;
    double *v;
    for (k=3;k<=14;++k){
        n = pow(2,k);
        v = generateV(n);
        sum = computeSum(v,n);
        printf("The difference in sums for N = %d is: %f \n",n, fabs(S - sum));
    }
    free(v);
    return 0;
}
