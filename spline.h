#ifndef III_SPLINE_H
#define III_SPLINE_H

#define PRECISION 4

typedef struct {
	int m, n;
} MatrixSize;

typedef struct {
	MatrixSize size;
	double **matrix;
} Matrix;

Matrix *add(Matrix *one, Matrix *two, double mult = 1, bool intoOne = false);

Matrix *mult(Matrix *one, Matrix *two, Matrix *res = nullptr);

double* spline(double data[][2], int pointNum);

Matrix *solveLU(Matrix *A, Matrix *b);

void LUfactor(Matrix *A, Matrix **L, Matrix **U, Matrix **P);

Matrix *substitute(Matrix *A, Matrix *b, bool front = true, Matrix *res = nullptr);

size_t matrixAllocSize(MatrixSize size);

Matrix* createDiagMatrix(MatrixSize size);

double** allocMatrix(MatrixSize size, bool zero = false);

Matrix* createMatrix(MatrixSize size, bool zero = false);

void destroyMatrix(Matrix *matrix);

Matrix* copyMatrix(Matrix* matrix);

bool cmpSize(Matrix *one, Matrix *two);

bool cmpSizeForMult(Matrix *one, Matrix *two);

void printMatrix(Matrix *matrix);

#endif //III_SPLINE_H
