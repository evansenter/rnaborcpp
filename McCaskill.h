// #define LA_COMPLEX_SUPPORT 1

#include <complex>

typedef std::complex<double> dcomplex;

std::complex<double>** runMcCaskill(char sequence[MAXSIZE], char structure[MAXSIZE]);
void solveZ(int i, int j, dcomplex x, char sequence[MAXSIZE], int *basePairs, int **basePairCounts, std::complex<double> **Z, std::complex<double> **ZB);
void solveZB(int i, int j, dcomplex x, char sequence[MAXSIZE], int *basePairs, int **basePairCounts, std::complex<double> **ZB, std::complex<double> **ZM);
void solveZM(int i, int j, dcomplex x, char sequence[MAXSIZE], int *basePairs, int **basePairCounts, std::complex<double> **ZB, std::complex<double> **ZM);
void solveLinearSystem(dcomplex **rootsOfUnity, dcomplex **Z);
int jPairedIn(int i, int j, int *basePairs);
int jPairedTo(int i, int j, int *basePairs);
int* getBasePairList(char *structure);
int numberOfBasePairs(int i, int j, int *basePairs);
int** fillBasePairCounts(int *basePairs, int n);
void printMatrix(dcomplex **matrix, char *title, int iStart, int iStop, int jStart, int jStop);
