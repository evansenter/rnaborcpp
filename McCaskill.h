double** runMcCaskill(char sequence[MAXSIZE]);
int solveZ( int i, int j, char sequence[MAXSIZE], char *structure, double **Z,  double **ZB, double *complex_x);
int solveZB(int i, int j, char sequence[MAXSIZE], char *structure, double **ZB, double **ZM, double *complex_x);
int solveZM(int i, int j, char sequence[MAXSIZE], char *structure, double **ZB, double **ZM, double *complex_x);
