#ifndef III_LAGRANGE_H
#define III_LAGRANGE_H

double* lagrange(double data[][2], int pointNum);

double calcDenominator(double data[][2], int pointNum, int index);

double calcMultiplier(double data[][2], int pointNum, int ignored, int level, int chosen = -1, double sum = 0, double product = 1);

#endif //III_LAGRANGE_H
