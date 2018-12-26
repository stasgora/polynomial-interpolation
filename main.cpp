#include <cstdio>
#include <fstream>
#include <limits>

#include "lagrange.h"
#include "spline.h"

using namespace std;

void pyplotPrint(double data[][2], int pointNum, double *lagParam, double *splParam, bool ignoreLagrange = false);

int main(int argc,  char* argv[]) {
	// dane, ilość punktów, co który na wejście
	if (argc < 4) {
		printf("Brak argumentów wywołania.");
		return 0;
	}
	int dataNum = stoi(argv[2]);
	int pointFreq = stoi(argv[3]);
	int pointNum = (dataNum + pointFreq - 1) / pointFreq;

	ofstream config;
	config.open ("../files/plot/plot_preset.py");
	config << "dataNum = " << dataNum << "\n";
	config << "freq = " << pointFreq << "\n";
	config.close();

	fstream dataStream;
	dataStream.open(argv[1], ios::in);
	if (!dataStream.is_open()) {
		printf("File open error");
		return 0;
	}
	double data[pointNum][2];
	long maxChars = numeric_limits<streamsize>::max();
	for (int i = 0; i < dataNum; ++i) {
		if(i % pointFreq == 0) {
			//data[i / pointFreq][0] = i;
			dataStream >> data[i / pointFreq][0];
			dataStream.ignore(maxChars, ',');
			dataStream >> data[i / pointFreq][1];
		}
		dataStream.ignore(maxChars, '\n');
	}
	double* lagParam;
	bool ignoreLagrange = false;
	if(!ignoreLagrange)
		lagParam = lagrange(data, pointNum);
	double* splParam = spline(data, pointNum);

	pyplotPrint(data, pointNum, lagParam, splParam, ignoreLagrange);

	free(lagParam);
	free(splParam);
	dataStream.close();
	return 0;
}

void pyplotPrint(double (*data)[2], int pointNum, double *lagParam, double *splParam, bool ignoreLagrange) {
	FILE *saved = stdout;
	if(!ignoreLagrange) {
		stdout = fopen("../files/plot/plot_lagrange.py", "w");

		printf("lagrange = [");
		for (int i = 0; i < pointNum; ++i) {
			printf("%.12e", lagParam[i]);
			if (i < pointNum - 1)
				printf(", ");
		}
		printf(" ]");

		fclose(stdout);
	}
	stdout = fopen("../files/plot/plot_spline.py", "w");

	printf("spline = [ ");
	for (int i = 0; i < pointNum - 1; ++i) {
		printf("[ ");
		for (int j = 0; j < 4; ++j) {
			printf("%.14f", splParam[i * 4 + j]);
			if (j < 3)
				printf(", ");
			else
				printf(" ]");
		}
		if (i < pointNum - 2)
			printf(", ");
	}
	printf(" ]");

	fclose(stdout);
	stdout = saved;
}
