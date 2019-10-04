#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "mkl_lapacke.h"

double *generate_matrix(int size)
{
    int i;
    double *matrix = (double *)malloc(sizeof(double) * size * size);
    srand(1);

    for (i = 0; i < size * size; i++)
    {
        matrix[i] = rand() % 100;
    }

    return matrix;
}

void print_matrix(const char *name, double *matrix, int size)
{
    int i, j;
    printf("matrix: %s \n", matrix);

    for (i = 0; i < size; i++)
    {
            for (j = 0; j < size; j++)
            {
                printf("%f ", matrix[i * size + j]);
            }
            printf("\n");
    }
}

int check_result(double *bref, double *b, int size) {
    int i;
    for(i=0;i<size*size;i++) {
        if (bref[i]!=b[i]) return 0;
    }
    return 1;
}
int my_dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb) {

    //Replace with your implementation
    LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv, b, ldb);
    
}


    void main(int argc, char *argv[])
    {

        int size = atoi(argv[1]);

        double *a, *aref;
        double *b, *bref;

        a = generate_matrix(size);
        aref = generate_matrix(size);        
        b = generate_matrix(size);
        bref = generate_matrix(size);
        

        //print_matrix("A", a, size);
        //print_matrix("B", b, size);

        // Using MKL to solve the system
        MKL_INT n = size, nrhs = size, lda = size, ldb = size, info;
        MKL_INT *ipiv = (MKL_INT *)malloc(sizeof(MKL_INT)*size);

        clock_t tStart = clock();
        info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
        printf("Time taken by MKL: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

        tStart = clock();    
        MKL_INT *ipiv2 = (MKL_INT *)malloc(sizeof(MKL_INT)*size);        
        my_dgesv(n, nrhs, a, lda, ipiv2, b, ldb);
        printf("Time taken by my implementation: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
        
        if (check_result(bref,b,size)==1)
            printf("Result is ok!\n");
        else    
            printf("Result is wrong!\n");
        
        //print_matrix("X", b, size);
        //print_matrix("Xref", bref, size);
    }
