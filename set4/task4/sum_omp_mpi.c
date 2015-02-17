#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "common.h"

Vector generateV(n){
    Vector v = createVector(n);//MPI(n, &WorldComm,1);
    double myNum;
    int i;
    for (i=1;i<v->len+1;++i){
        myNum = 1/(i*i*1.0);
        v->data[i-1] = myNum;
    }
    return v;
}


int main(int argc, char **argv)
{ 
  // Define variables
  int n, rank, size, i;
  double S = M_PI*M_PI/6;
  double sum, t1, t2, dt, error;
  Vector v;
  int *len, *displ;
  double *myData;
    
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  int k;
  for (k = 3; k<=14; ++k){
  n = pow(2,k);
  int elements_per_proc = n/size;
  if (rank == 0) {    
    //n = elements_per_proc*size;
    v = generateV(n);
    myData = v->data;
    t1 = MPI_Wtime();
  }
  
  // buffer for data
  double *privateData = malloc(sizeof(double) * elements_per_proc);
  MPI_Scatter(myData, elements_per_proc, MPI_DOUBLE, privateData,
            elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Compute the sum
  sum = 0.0;
  #pragma omp parallel for schedule(static) reduction(+:sum)  
  for (i = 0; i < elements_per_proc; ++i) {
     sum = sum + privateData[i];
  }
  
  double *ourSum = NULL;
  if (rank == 0) {
    ourSum = malloc(sizeof(double) * size);
  }
  
  MPI_Gather(&sum, 1, MPI_DOUBLE, ourSum, 1, MPI_DOUBLE, 0,
           MPI_COMM_WORLD);


  if (rank == 0) {
    double allSums = 0;
    for (i = 0; i < size; ++i) {
      allSums = allSums + ourSum[i];
    }
    t2 = MPI_Wtime();
    dt = t2 - t1;
    error = fabs(S-allSums);
    printf ("sum= %f error= %f dt= %f n=%d \n", allSums, error, dt, n);
  }
  }
  MPI_Finalize();
  return 0;
}
