#include <cstdlib>
#include <cstdio>
#include "lagrange.h"

double *lagrange(double data[][2], int pointNum) {
	double *params = (double*) malloc(sizeof(double) * pointNum);
	for (int i = 0; i < pointNum; ++i)
		params[i] = 0;
	for (int i = 0; i < pointNum; ++i) {
		double denominator = calcDenominator(data, pointNum, i);
		//printf("%f ", denominator);
		for (int j = 0; j < pointNum; ++j) {
			params[j] += (data[i][1] * calcMultiplier(data, pointNum, i, j)) * (j % 2 == 0 ? 1 : -1) / denominator;
			//printf("%f ", data[i][1] * calcMultiplier(data, pointNum, i, j));
		}
		//printf("\n");
	}
	return params;
}

double calcDenominator(double data[][2], int pointNum, int index) {
	double res = 1;
	for (int i = 0; i < pointNum; ++i)
		if (i != index)
			res *= data[index][0] - data[i][0];
	return res;
}

double calcMultiplier(double data[][2], int pointNum, int ignored, int level, int chosen, double sum, double product) {
	chosen >= 0 ? product *= data[chosen][0] : level++;
	if(level <= 1)
		return sum + product;
	for (int i = chosen + 1; i < pointNum - level + 2; ++i)
		if (i != ignored)
			sum = calcMultiplier(data, pointNum, ignored, level - 1, i, sum, product);
	return sum;
}
