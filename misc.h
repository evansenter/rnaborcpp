int ElementInList(char subshape[SHAPELENGTH],char shapeList[SHAPELENGTH][SHAPELENGTH],int listlen);
int NumInList(int a, int L[SHAPELENGTH],int listlen);
double MaxInArray(double array[],int arraylen);
int BP(int i,int j,char sequence[MAXSIZE]);
double ***Allocate3DMatrix(int a, int b, int c);
double **Allocate2DMatrix(int a, int b);
int CheckSequence(char sequence[MAXSIZE]);
int Free3DMatrix(double ***Matrix, int a, int b, int c);
int Free2DMatrix(double **Matrix,int a, int b);
