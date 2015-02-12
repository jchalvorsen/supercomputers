#include <stdio.h>
#include <math.h>
#include "common.h"


Vector generateV(n){
    Vector v = createVector(n);
    double myNum;
    int i;
    for (i=1;i<v->len+1;++i){
        myNum = 1/(i*i*1.0);
        v->data[i-1] = myNum;
    }
    return v;    
}

void printVector(Vector v){
    int i;
    for (i=0;i<v->len;++i){
        printf("%f \n", v->data[i]);
    }
}

double computeSum(Vector v){
    double s = 0;
    int i;
    for (i=0;i<v->len;++i){
        s += v->data[i];
    }
    return s;
}

int main(int argc, char **argv )
{    
    // Loop over and print the sum differences:
    int k, n;
    double S = M_PI*M_PI/6;
    double sum;
    Vector v;
    for (k=3;k<=14;++k){
        n = pow(2,k);
        v = generateV(n);
        if (k==3)
            printVector(v);
        sum = computeSum(v);
        printf("The difference in sums for N = %d is: %f \n",n, fabs(S - sum));
    }
    
    return 0;
}



