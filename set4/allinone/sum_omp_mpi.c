#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
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

int main(int argc, char **argv)
{ 
  // Define variables
  int k, n, rank, size, i;
  double S = M_PI*M_PI/6, sum = 0.0;
  double t1, t2, dt, error;
  double *v = NULL;
  int *len, *displ;
  
  // initialize MPI and get arguments  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  for (k = 3; k<=14; ++k){
  n = pow(2,k);
  int elements_per_proc = n/size;
  if (rank == 0) {    
    v = generateV(n);
    t1 = MPI_Wtime();
  }
  
  // buffer for data
  double *privateData = malloc(sizeof(double) * elements_per_proc);
  MPI_Scatter(v, elements_per_proc, MPI_DOUBLE, privateData,
            elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Compute the sum (on all processors)
  #pragma omp parallel for schedule(static) reduction(+:sum)
  for (i = 0; i < elements_per_proc; ++i) {
     sum = sum + privateData[i];
  } 
  
  // Gather the result, using the reduce function.
  double *globalsum = malloc(sizeof(double));
  MPI_Reduce(&sum,globalsum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  if (rank == 0) {
    t2 = MPI_Wtime();
    dt = t2 - t1;
    error = fabs(S-globalsum[0]);
    printf ("sum= %f error= %f dt= %f n=%d \n", globalsum[0], error, dt, n);
  }
  }
  MPI_Finalize();
  return 0;
}
