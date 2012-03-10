#define LA_COMPLEX_SUPPORT 1

#include <complex>

typedef std::complex<double> dcomplex;

void solveZB(int i, int j, dcomplex x, char sequence[MAXSIZE], int **basePairCounts, std::complex<double> **ZB, std::complex<double> **ZM);
void solveZM(int i, int j, dcomplex x, char sequence[MAXSIZE], int **basePairCounts, std::complex<double> **ZB, std::complex<double> **ZM);
void solveZ(int i, int j, dcomplex x, char sequence[MAXSIZE], int **basePairCounts, std::complex<double> **Z, std::complex<double> **ZB);
std::complex<double>** runMcCaskill(char sequence[MAXSIZE]);
int* getBasePairList(char *structure);
int numberOfBasePairs(int i, int j, int *basePairs);
int** fillBasePairCounts(int *basePairs, int n);