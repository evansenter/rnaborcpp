double** runMcCaskill(char *sequence);
int solveZ(int i, int j, int x, char *sequence, char *structure, double ***Z, double ***ZB, double **complexRoots);
int solveZB(int i, int j, int x, char *sequence, char *structure, double ***ZB, double ***ZM, double **complexRoots);
int solveZM(int i, int j, int x, char *sequence, char *structure, double ***ZB, double ***ZM, double **complexRoots);
double* addComplex(double *complexA, double *complexB, double *complexSum);
double* scalarComplex(double *complex, double scalar, double *complexScalar);
double* productComplex(double *complexA, double *complexB, double *complexProduct);
double* powerComplex(double *complex, int power, double *complexPower);
double* rectangularRootComplex(double angle, double *complexResult);