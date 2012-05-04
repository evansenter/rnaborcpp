double** runMcCaskill(char sequence[MAXSIZE]);
void solveZ(int i, int j, dcomplex x, char sequence[MAXSIZE], int *basePairs, int **basePairCounts, double **Z);
void solveLinearSystem(double **rootsOfUnity);
int jPairedIn(int i, int j, int *basePairs);
int jPairedTo(int i, int j, int *basePairs);
int* getBasePairList(char *structure);
int numberOfBasePairs(int i, int j, int *basePairs);
int** fillBasePairCounts(int *basePairs, int n);
void printMatrix(double **matrix, char *title, int iStart, int iStop, int jStart, int jStop);