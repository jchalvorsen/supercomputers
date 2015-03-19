/*
  C-program to solve the two-dimensional Poisson equation on
  a unit square using one-dimensional eigenvalue decompositions
  and fast sine transforms
  einar m. ronquist
  ntnu, october 2000
  revised, october 2001
*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

typedef double Real;

/* function prototypes */
Real *createRealArray (int n);
Real **createReal2DArray (int m, int n);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);
void printMatrix(Real **m, int x, int y);
void printVector(Real *v, Real n);
void transpose(Real **b, int *numberOfCols, int *startCol, Real *sendbuffer, Real *rbuffer);

int main(int argc, char **argv )
{
    Real *diag, **b, *z, *sendbuffer, *rbuffer;
    Real pi, h, umax;
    int i, j, n, m, nn , rank, size, rest, *numberOfCols,*startCol,tcount;

    // initialize MPI and get arguments
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* the total number of grid points in each spatial direction is (n+1) */
    /* the total number of degrees-of-freedom in each spatial direction is (n-1) */
    /* this version requires n to be a power of 2 */
    if( argc < 2 ) {
        printf("need a problem size\n");
        return 1;
    }

    n  = atoi(argv[1]);
    m  = n-1;


    // Deciding what processor gets what data.
    numberOfCols = malloc( size * sizeof(int) );
    rest = m % size;
    for (i = 0; i < size-rest; ++i){
        numberOfCols[i] = m/size;
    }
    for (i = size-rest; i < size; ++i){
        numberOfCols[i] = m/size+1;
    }
    startCol = malloc( size * sizeof(int) );
    startCol[0] = 0;

    for (i = 1; i < size; ++i){
        startCol[i] = startCol[i-1] + numberOfCols[i-1];
    }

    //                      //
    // Assigning constants  //
    //                      //
    h    = 1./(Real)n;
    pi   = 4.*atan(1.);

    // Buffer for sine transform
    nn = 4*n;
    z    = createRealArray (nn);

    // buffers for transpose
    sendbuffer = createRealArray (m*numberOfCols[rank]);
    rbuffer = createRealArray (m*numberOfCols[rank]);

    // b will have a different sizes given thread
    b    = createReal2DArray (numberOfCols[rank],m);

    // diag should be available to all threads
    diag = createRealArray (m);
    #pragma omp parallel for
    for (i=0; i < m; i++) {
        diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
    }

    // Fill in values in b. TODO: map from loading function
    #pragma omp parallel for private(i)
    for (j=0; j < numberOfCols[rank]; j++) {
        for (i=0; i < m; i++) {
            b[j][i] = h*h;//i + m*j + startCol[rank]*m + 1;
        }
    }

    //                      //
    // Start of algorithm   //
    //                      //
     
    
    #pragma omp parallel shared(b) private(z)
    {
	    z    = createRealArray (nn);
	    #pragma omp for 
	    for (i=0; i < numberOfCols[rank]; i++) {
		fst_(b[i], &n, z, &nn);
	    }
    }
    transpose(b, numberOfCols, startCol, sendbuffer,rbuffer);

    #pragma omp parallel shared(b) private(z)
    {
	    z = createRealArray (nn);
  	    #pragma omp for 
	    for (i=0; i < numberOfCols[rank]; i++) {
		fstinv_(b[i], &n, z, &nn);
	    }
    }

    // Find the eigenvalues (b is now transposed, but does not matter here)
    #pragma omp parallel for private(i)
    for (j=0; j < numberOfCols[rank]; j++) {
        for (i=0; i < m; i++) {
            b[j][i] = b[j][i]/(diag[i] + diag[j + startCol[rank]]);
        }
    }
    

    //Comment: Do not use all threads..!
    #pragma omp parallel shared(b) private(z)
    {
	    z    = createRealArray (nn);
	    #pragma omp for 
	    for (i=0; i < numberOfCols[rank]; i++) {
		fst_(b[i], &n, z, &nn);
	    }
    }

    transpose(b, numberOfCols, startCol, sendbuffer,rbuffer);

    #pragma omp parallel shared(b) private(z)
    {
	    z    = createRealArray (nn);
	    #pragma omp for
	    for (i=0; i < numberOfCols[rank]; i++) {
		fstinv_(b[i], &n, z, &nn);
	    }
    }  
    //                      //
    // End of algorithmnts  //
    //                      //

    // Get max and max-reduce
    umax = 0.0;
    for (j=0; j < numberOfCols[rank]; j++) {
        for (i=0; i < m; i++) {
            if (b[j][i] > umax) umax = b[j][i];
        }
    }
    double *globalsum = malloc(sizeof(double));
    MPI_Reduce(&umax,globalsum,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if(rank == 0){
        printf (" umax = %e \n",globalsum[0]);
        // Print to file:
        FILE *f = fopen("out.txt", "a");
        if (f == NULL)
        {
            printf("Error opening file!\n");
            exit(1);
        }

        /* print some text */
        fprintf(f,"umax = %e \n",globalsum[0]);

        fclose(f);
    }
    // MPI IO print to file:

//    MPI_File fh;
//    MPI_File_open(MPI_COMM_WORLD, "out2.dat", MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

//    MPI_Offset offset = startCol[rank]*m;

//    MPI_File_iwrite_at(fh, offset, b[0], numberOfCols[rank]*m, MPI_DOUBLE, MPI_STATUS_IGNORE);
//    MPI_File_close(&fh);

    free(diag);
    free(b);
    free(z);
    free(sendbuffer);
    free(rbuffer);
    free(globalsum);
    free(numberOfCols);
    free(startCol);
    MPI_Finalize();
    return 0;
}

void transpose(Real **b, int *numberOfCols, int *startCol, Real *sendbuffer, Real *rbuffer)
{
    // Transpose the matrix b in place over MPI

    // Lay out all elements in right memory order
    int count = 0;
    int p, i , j, rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for (p = 0; p < size; ++p){
        for (i = 0; i < numberOfCols[rank]; ++i){
            for (j = 0; j < numberOfCols[p]; ++j){
                sendbuffer[count] = b[i][j + startCol[p]];
                count += 1;
            }
        }
    }

    int *srdispls = malloc( size * sizeof(int) );
    int *srcounts = malloc( size * sizeof(int) );
    for (i=0; i < size; ++i){
        srcounts[i] = numberOfCols[i]*numberOfCols[rank];
        srdispls[i] = startCol[i]*numberOfCols[rank];
    }
    MPI_Alltoallv(sendbuffer, srcounts, srdispls, MPI_DOUBLE, rbuffer, srcounts, srdispls, MPI_DOUBLE, MPI_COMM_WORLD);

    // Taking the data back to b
    count = 0;
    for (p = 0; p < size; ++p){
        for (i = 0; i < numberOfCols[p]; ++i){
            for (j = 0; j < numberOfCols[rank]; ++j){
                b[j][i + startCol[p]] = rbuffer[count];
                count += 1;
            }
        }
    }
}

void printMatrix(Real **m, int x, int y){
    int i, j;
    for (j=0; j < y; j++) {
        for (i=0; i < x; i++) {
            printf("|%2.2f ",m[i][j]);
        }
        printf("| \n");
    }
}

void printVector(Real *v, Real n){
    int i;
    for (i = 0; i < n; i++){
        printf("%d \n", v[i]);
    }
}

Real *createRealArray (int n)
{
    Real *a;
    int i;
    a = (Real *)malloc(n*sizeof(Real));
    for (i=0; i < n; i++) {
        a[i] = 0.0;
    }
    return (a);
}

Real **createReal2DArray (int n1, int n2)
{
    int i, n;
    Real **a;
    a    = (Real **)malloc(n1   *sizeof(Real *));
    a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
    for (i=1; i < n1; i++) {
        a[i] = a[i-1] + n2;
    }
    n = n1*n2;
    memset(a[0],0,n*sizeof(Real));
    return (a);
}
