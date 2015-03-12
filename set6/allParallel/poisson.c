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
void transpose (Real **bt, Real **b, int m);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);
void procedure (Real **bt, Real **b, int n, int m, Real *z, int nn);
void printMatrix(Real **m, int x, int y);
void printVector(Real *v, int n);


int main(int argc, char **argv )
{
    Real *diag, **b, **bt, *z;
    Real pi, h, umax;
    int i, j, n, m, nn, rank, size;

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
    int *numberOfCols = malloc( size * sizeof(int) );
    numberOfCols[0] = m/size;
    if (size > 1){
        for (i = 1; i < size; ++i){
            numberOfCols[i] = n/size;
        }
    }
    int *startCol = malloc( size * sizeof(int) );
    startCol[0] = 0;

    if (size > 1){
        startCol[1] = m/size;
        if (size > 2){
            for (i = 2; i < size; ++i){
                startCol[i] = startCol[i-1] + n/size;
            }
        }
    }

    // Assigning constants
    // b will have a different sizes given thread
    b    = createReal2DArray (numberOfCols[rank],m);

    h    = 1./(Real)n;
    pi   = 4.*atan(1.);

    // diag should be available to all threads
    diag = createRealArray (m);
    for (i=0; i < m; i++) {
        diag[i] = 2.*(1.-cos((i+1)*pi/(Real)n));
    }

    // Fill in values in b. TODO: map from loading function
    for (j=0; j < numberOfCols[rank]; j++) {
        for (i=0; i < m; i++) {
            b[j][i] = i + m*j + startCol[rank]*m + 1;
        }
    }
    printf("Proc number %d says hi: \n", rank);
    printMatrix(b, numberOfCols[rank], m);

    MPI_Finalize();
    return 0;
}

void procedure (Real **bt, Real **b, int n, int m, Real *z, int nn){
    int i, size, rank;
    Real **test;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //printf("start of procedure \n");

    // All processors need this internal memory managed
    test = createReal2DArray(m/size,m);
    int *srdispls = malloc( size * sizeof(int) );
    int *srcounts = malloc( size * sizeof(int) );
    for (i=0; i < size; ++i){
        srcounts[i] = m*m/size;
        srdispls[i] = i*m*m/size;
    }
    if (rank == 0){

        for (i=0; i < size; ++i){
            printf("%f \n",b[0][srdispls[i]]);
        }
        printf("This is whole matrix \n");
        printMatrix(b, m, m);
        printf("%p \n", &(b[0][0]));
        printf("%p \n", (b[0]));
        printf("This is whole vector \n");
        printVector(b[0], m*m);
    }

    //printf("scattering begins on p = %d \n",rank);
    MPI_Scatterv(&b[0], srcounts, srdispls,
            MPI_DOUBLE, &test[0], m*m/size,
            MPI_DOUBLE,
            0, MPI_COMM_WORLD);
    //printf("This is proc nr %d's matrix \n", rank);
    //printVector(test[0], m/size*m);

    // print local data:
    int p;
    for (p=0; p<size; p++) {
        if (rank == p) {
            printf("Local process on rank %d is:\n", rank);
            printMatrix(test, m/size, m);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }


    // First do the transform on all the coloumns
    for (i=0; i < m/size; i++) {
        fst_(test[i], &n, z, &nn);
    }
    printf("This is proc nr %d's matrix after fst %f \n", rank, test[0][0]);
    printMatrix(test, m/size, m);
    MPI_Gatherv(test[0], m*m/size, MPI_DOUBLE,
            b[0], srcounts, srdispls,
            MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0){
        printf("This is whole reassambled matrix \n");
        printMatrix(b, m, m);
    }
    if (rank == 0){
        // Transpose it
        transpose (bt,b,m);

        // Then do the inverse transform.
        for (i=0; i < m; i++) {
            fstinv_(bt[i], &n, z, &nn);
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

void printVector(Real *v, int n){
    int i;
    for (i = 0; i < n; i++){
        printf("%2.2f \n", v[i]);
    }
}


void transpose (Real **bt, Real **b, int m)
{
    int i, j;
    for (j=0; j < m; j++) {
        for (i=0; i < m; i++) {
            bt[j][i] = b[i][j];
        }
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