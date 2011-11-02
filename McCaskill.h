double** runMcCaskill(char *sequence);
int solveZ(int i, int j, int x, char *sequence, char *structure, double ***Z, double ***ZB, double **complexRoots);
int solveZB(int i, int j, int x, char *sequence, char *structure, double ***ZB, double ***ZM, double **complexRoots);
int solveZM(int i, int j, int x, char *sequence, char *structure, double ***ZB, double ***ZM, double **complexRoots);
void writeScalarProduct(double ***matrix, int i, int j, int realIndex, int complexIndex, double scalar);