#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "mkl_lapacke.h"

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
    printf("\n\nmatrix %s : \n \n", name);

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
    //print_matrix("bref", bref, size);
    //print_matrix("b", b, size);
    for(i=0;i<size*size;i++) {
        if (abs(bref[i] - b[i]) > 0.01 )
        {
        	printf("The first difference found is at the %ith iteration : \n",i+1);
        	printf("b_%i = %5.30f ",i,b[i]);
        	printf("bref_%i = %5.30f ",i,bref[i]);
        	return 0;	
        } 
    }
    return 1;
}

double *matrixProduct(double *leftMatrix, double *rightMatrix, int size)
{
	double *XMatrix = (double *)malloc(sizeof(double) * size * size);

	for(int i = 0; i < size * size; ++i)
			XMatrix[i] = 0.0;

	for(int i = 0; i < size; ++i)
		for(int j = 0; j < size; ++j)
			for (int k = 0; k < size; ++k)
				XMatrix[i * size + j] += leftMatrix[i * size + k] * rightMatrix[k * size + j];
		
	return XMatrix;
}

void QRdecomposition(double *matrix, double *QMatrix, double *RMatrix, int size)
{
	for (int i = 0; i < size * size; ++i)
	{
    	RMatrix[i] = 0.0;
		QMatrix[i] = 0.0;
	}
	double percentage = 0.0;

	for(int i = 0; i < size; ++i)
	{
		printf("\r");
		printf("QRdecomposition process: %u %%", (int) (percentage*100));
		percentage += 1.0 / size;

		double RMatrix_i_j = 0.0;
		
		for(int j = 0; j < size; ++j)
			RMatrix_i_j += matrix[j * size + i] * matrix[j * size + i];

		RMatrix[i * (size + 1)] = sqrt(RMatrix_i_j);
		
		for(int j = 0; j < size; ++j)
			QMatrix[j * size + i] = matrix[j * size + i] / RMatrix[i * (size + 1)];
		
		for(int j = i+1; j < size; ++j)
		{
			RMatrix_i_j = 0;

			for(int k = 0; k < size; ++k)
				RMatrix_i_j += matrix[k*size + j] * QMatrix[k*size + i];
			
			RMatrix[i*size + j] = RMatrix_i_j;

			for(int k = 0; k < size; ++k)
				matrix[k * size + j] -= RMatrix[i * size + j] * QMatrix[k * size + i];
		}
	}
	printf("\r");
	printf("QRdecomposition process: 100 %%");
}

void solveQRXeqB(double *QMatrix, double *RMatrix, double *XMatrix, double *BMatrix, int size)
{
	for (int i = 0; i < size * size; ++i)
    	XMatrix[i] = 0.0;

    double *transposedQMatrix = (double *)malloc(sizeof(double) * size * size);

	for(int i = 0; i < size; ++i) 
   		for(int j = 0; j < size; ++j) 
      		transposedQMatrix[i * size + j] = QMatrix[j * size + i];

	double *transposedQBMatrix = matrixProduct(transposedQMatrix, BMatrix, size);

	double percentage = 0.0;

	printf("\n");

	for(int i = 0; i < size; ++i) {
		printf("\r");
		printf("Solving: %u %%", (int)(percentage * 100));
		percentage += 1.0 / size;

		for(int j = 0; j < size; ++j) {
			double sou = 0.0;

			for(int sub = size - i; sub < size; ++sub){
				sou += RMatrix[(size - i - 1) * size + sub]*XMatrix[sub * size + (size - j - 1)];
			}
			XMatrix[(size - i - 1) * size + (size - j - 1)] = (transposedQBMatrix[(size - i - 1) * size + (size - j - 1)] - sou)/RMatrix[(size - i - 1) * (size+1)];
		}
	}
	printf("\r");
	printf("Solving: 100 %%");
}

void main(int argc, char *argv[])
{

    int size = atoi(argv[1]);
    //int size = 500;

    double *a, *aref;
    double *b, *bref;

    a = generate_matrix(size);
    aref = generate_matrix(size);        
    b = generate_matrix(size);
    bref = generate_matrix(size);

    //print_matrix("A", a, size);
    //print_matrix("B", b, size);

    // Using MKL to solve the system
    int n = size, nrhs = size, lda = size, ldb = size, info;
    //int *ipiv = (MKL_INT *)malloc(sizeof(MKL_INT)*size);

	clock_t tStart;

    //info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
    //printf("Time taken by MKL: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

    tStart = clock();



    double *Q = (double *)malloc(sizeof(double) * size * size);
    double *R = (double *)malloc(sizeof(double) * size * size);

    QRdecomposition(a, Q, R, size);



 	double *x = (double *)malloc(sizeof(double) * size * size);

	solveQRXeqB(Q, R, x, b, size);
	b = matrixProduct(matrixProduct(Q,R,size),x,size);
    //print_matrix("X", x, size);
    //print_matrix("B", b, size);


    printf("\nTime : %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
    
    if (check_result(bref,b,size)==1)
        printf("Result is ok!\n");
    else    
        printf("Result is wrong!\n");
}
