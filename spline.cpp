#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "spline.h"

#define SPLINE_PARAM_NUM(X) (4 * ((X) - 1))

double *spline(double data[][2], int pointNum) {
	int paramNum = SPLINE_PARAM_NUM(pointNum);
	Matrix *b = createMatrix({paramNum, 1}, true);
	Matrix *A = createMatrix({paramNum, paramNum}, true);

	double x_diff = data[1][0] - data[0][0];

	for (int i = 0; i < pointNum - 1; ++i) {
		A->matrix[i * 2][i * 4] = 1;
		b->matrix[i * 2][0] = data[i][1];

		double x_val = 1;
		for (int j = 0; j < 4; ++j) {
			A->matrix[i * 2 + 1][i * 4 + j] = x_val;
			x_val *= x_diff;
		}
		b->matrix[i * 2 + 1][0] = data[i + 1][1];
	}
	int m_off = 2 * pointNum - 2;
	for (int i = 0; i < pointNum - 2; ++i) {
		A->matrix[m_off + i * 2][i * 4 + 1] = 1; //b
		A->matrix[m_off + i * 2][i * 4 + 2] = 2 * x_diff; //c
		A->matrix[m_off + i * 2][i * 4 + 3] = 3 * x_diff * x_diff; //d
		A->matrix[m_off + i * 2][i * 4 + 5] = -1; //b+1

		A->matrix[m_off + i * 2 + 1][i * 4 + 2] = 2; //c
		A->matrix[m_off + i * 2 + 1][i * 4 + 3] = 6 * x_diff; //d
		A->matrix[m_off + i * 2 + 1][i * 4 + 6] = -2; //c+1
	}
	m_off += 2 * pointNum - 4;
	A->matrix[m_off    ][2] = 1; //c
	A->matrix[m_off + 1][(pointNum - 2) * 4 + 2] = 2; //c
	A->matrix[m_off + 1][(pointNum - 2) * 4 + 3] = 6 * x_diff; //d

	Matrix *x = solveLU(A, b);
	double *params = (double*) malloc(sizeof(double) * paramNum);
	for (int i = 0; i < paramNum; ++i)
		params[i] = x->matrix[i][0];
	destroyMatrix(x);
	destroyMatrix(A);
	destroyMatrix(b);
	return params;
}

Matrix *solveLU(Matrix *A, Matrix *b) {
	Matrix *L, *U, *P;
	LUfactor(A, &L, &U, &P);

	Matrix *b_2 = mult(P, b);
	Matrix *y = substitute(L, b_2);
	Matrix *x = substitute(U, y, false);
	destroyMatrix(L);
	destroyMatrix(U);
	destroyMatrix(P);
	destroyMatrix(y);
	destroyMatrix(b_2);

	return x;
}

void LUfactor(Matrix *A, Matrix **L, Matrix **U, Matrix **P) {
	*U = copyMatrix(A);
	*L = createDiagMatrix(A->size);
	*P = createDiagMatrix(A->size);

	for (int k = 0; k < A->size.m - 1; ++k) {
		int maxIndex = k;
		double max = fabs((*U)->matrix[k][k]);
		for (int i = k + 1; i < A->size.m; ++i)
			if(fabs((*U)->matrix[i][k]) > max) {
				max = fabs((*U)->matrix[i][k]);
				maxIndex = i;
			}
		double temp;
		for (int i = k; i < A->size.m; ++i) {
			temp = (*U)->matrix[k][i];
			(*U)->matrix[k][i] = (*U)->matrix[maxIndex][i];
			(*U)->matrix[maxIndex][i] = temp;
		}
		for (int i = 0; i < k; ++i) {
			temp = (*L)->matrix[k][i];
			(*L)->matrix[k][i] = (*L)->matrix[maxIndex][i];
			(*L)->matrix[maxIndex][i] = temp;
		}
		double *row = (*P)->matrix[k];
		(*P)->matrix[k] = (*P)->matrix[maxIndex];
		(*P)->matrix[maxIndex] = row;

		for (int j = k + 1; j < A->size.m; ++j) {
			(*L)->matrix[j][k] = (*U)->matrix[j][k] / (*U)->matrix[k][k];
			for (int i = k; i < A->size.m; ++i)
				(*U)->matrix[j][i] -= (*L)->matrix[j][k] * (*U)->matrix[k][i];
		}
	}
}

Matrix* add(Matrix *one, Matrix *two, double mult, bool intoOne) {
	if (!cmpSize(one, two)) {
		printf("matrix sizes don't match\n");
		return nullptr;
	}
	Matrix *res = nullptr;
	if (!intoOne)
		res = createMatrix(one->size);
	for (int i = 0; i < one->size.m; ++i)
		for (int j = 0; j < one->size.n; ++j)
			if (!intoOne)
				res->matrix[i][j] = one->matrix[i][j] + mult * two->matrix[i][j];
			else
				one->matrix[i][j] = one->matrix[i][j] + mult * two->matrix[i][j];
	if (!intoOne)
		return res;
	else
		return one;
}

Matrix *mult(Matrix *one, Matrix *two, Matrix *res) {
	if (!cmpSizeForMult(one, two) || (res != nullptr && (one->size.m != res->size.m || two->size.n != res->size.n))) {
		printf("matrix sizes don't match\n");
		return nullptr;
	}
	if (res == nullptr)
		res = createMatrix({one->size.m, two->size.n});
	for (int i = 0; i < one->size.m; ++i)
		for (int j = 0; j < two->size.n; ++j) {
			double sum = 0;
			for (int k = 0; k < one->size.n; ++k)
				sum += one->matrix[i][k] * two->matrix[k][j];
			res->matrix[i][j] = sum;
		}
	return res;
}

Matrix *substitute(Matrix *A, Matrix *b, bool front, Matrix *res) {
	int start = front ? 0 : A->size.m - 1, end = front ? A->size.m : 0;
	if (res == nullptr)
		res = createMatrix(b->size);
	for (int i = start; front ? i < end : i >= end; i += front ? 1 : -1) {
		double sum = b->matrix[i][0];
		int start2 = front ? 0 : i + 1, end2 = front ? i : A->size.m;
		for (int j = start2; j < end2; j++)
			sum -= res->matrix[j][0] * A->matrix[i][j];
		res->matrix[i][0] = sum / A->matrix[i][i];
	}
	return res;
}

void allocMatrix(Matrix *matrix, bool zero) {
	if (zero)
		matrix->matrix = (double **) calloc(matrixAllocSize(matrix->size), 1);
	else
		matrix->matrix = (double **) malloc(matrixAllocSize(matrix->size));
	if (matrix->matrix == nullptr) {
		printf("failed to alloc matrix");
		return;
	}
	for (int i = 0; i < matrix->size.m; ++i) {
		matrix->matrix[i] = (double *) matrix->matrix + matrix->size.m + i * matrix->size.n;
	}
}

size_t matrixAllocSize(MatrixSize size) {
	return sizeof(double) * size.m * size.n + sizeof(double *) * size.m;
}

Matrix *createMatrix(MatrixSize size, bool zero) {
	Matrix * matrix = (Matrix*) malloc(sizeof(Matrix));
	matrix->size = size;
	allocMatrix(matrix, zero);
	return matrix;
}

Matrix* createDiagMatrix(MatrixSize size) {
	Matrix *diag = createMatrix(size, true);
	for (int i = 0; i < size.m; ++i)
		diag->matrix[i][i] = 1;
	return diag;
}

void destroyMatrix(Matrix *matrix) {
	free(matrix->matrix);
	free(matrix);
}

Matrix* copyMatrix(Matrix *matrix) {
	Matrix *copy = createMatrix(matrix->size);
	for (int i = 0; i < matrix->size.m; ++i)
		for (int j = 0; j < matrix->size.n; ++j)
			copy->matrix[i][j] = matrix->matrix[i][j];
	return copy;
}

bool cmpSize(Matrix *one, Matrix *two) {
	return one->size.m == two->size.m && one->size.n == two->size.n;
}

bool cmpSizeForMult(Matrix *one, Matrix *two) {
	return one->size.n == two->size.m;
}

void printMatrix(Matrix *matrix) {
	for (int i = 0; i < matrix->size.m; ++i) {
		printf("| ");
		for (int j = 0; j < matrix->size.n; ++j) {
			printf(" %.*f", matrix->matrix[i][j] < 0 ? PRECISION - 1 : PRECISION, matrix->matrix[i][j]);
		}
		printf(" |\n");
	}
	printf("\n");
}
