// clang++ rnabor.cpp -I /usr/local/include/lapackpp `pkg-config lapackpp --libs` -o rnabor

#include <lapackpp.h>

#define LA_COMPLEX_SUPPORT 1
#define NBPAIRS 7
#define MAXLOOP 30
#define MAXSIZE 1100
#define MAXALPHA 20 /* maximal length of alphabet */
#define MIN_PAIR_DIST 3
#define MAX_INTERIOR_DIST 30
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define MIN2(A, B) ((A) < (B) ? (A) : (B))
#define SET_Z(i, j, value) \
Z[i][j] = value; \
Z[j][i] = value;

int seqlen;
int energy_set = 0; /* 0 = BP; 1 = any with GC; 2 = any with AU parameters */
short *S, *S0, *S1;
short alias[MAXALPHA+1];
static const char Law_and_Order[] = "_ACGUTXKI";
int pair[MAXALPHA+1][MAXALPHA+1];
int tetra_loop;
int MAX_NINIO = 300;
double kT;
double ML_base;
double ML_close;
double temperature = 37.0;

double **Allocate2DMatrix(int m, int n);
int canBasePair(int i, int j, char sequence[MAXSIZE]);
double** runMcCaskill(char sequence[MAXSIZE]);
void solveZ(int i, int j, char sequence[MAXSIZE], double **Z, double **ZB);
void solveZB(int i, int j, char sequence[MAXSIZE], double **ZB, double **ZM);
void solveZM(int i, int j, char sequence[MAXSIZE], double **ZB, double **ZM);

double **Allocate2DMatrix(int m, int n) {
  int i;
  double **Matrix;
  Matrix = (double **) calloc(m, sizeof (double *));
  
  if (Matrix == NULL) {
    printf("Out of memory.\n");
    exit(1);
  }
  
  for (i = 0; i < m; ++i) {
    Matrix[i] = (double *) calloc(n, sizeof(double));
    
    if (Matrix[i] == NULL) {
      printf("Out of memory.\n");
      exit(1);
    }
  }
  
  return Matrix;
}

int canBasePair(int i, int j, char sequence[MAXSIZE]) {
  if (j - i <= MIN_PAIR_DIST) {
    return 0;
  } else if ((sequence[i] == 'A' && sequence[j] == 'U') || (sequence[i] == 'U' && sequence[j] == 'A')) {
    return 1;
  } else if ((sequence[i] == 'C' && sequence[j] == 'G') || (sequence[i] == 'G' && sequence[j] == 'C')) {
    return 1;
  } else if ((sequence[i] == 'U' && sequence[j] == 'G') || (sequence[i] == 'G' && sequence[j] == 'U')) {
    return 1;
  } else {
    return 0;
  }
}

int CheckSequence(char sequence[MAXSIZE]){     //safe
  int i;
  for (i=1;i<strlen(sequence);++i)
    {
     if (toupper(sequence[i])!='A' && toupper(sequence[i])!='U' && toupper(sequence[i])!='C' && toupper(sequence[i])!='G'&& toupper(sequence[i])!='T') //check if there are invalid characters
        {
         printf("The input string should only contain A,U,T,C,G!\n");
         exit(1);
        }
     else if (toupper(sequence[i])=='T') //change T to U
        sequence[i]='U';
     else
        sequence[i]=toupper(sequence[i]); // change lower case to upper case
    }
  return 0;
}

void Initialize_Params() {
  P = scale_parameters();//from params.c, gets our parameters, for a given temp
}

int main(int argc, char *argv[]) {
  char sequence[] = " GGGGGCCCCCGGGGGCCCCCGGGGGCCCCC"; // Sequences start at index 1
  CheckSequence(sequence);
  S0=encode_seq(sequence + 1);
  seqlen=strlen(sequence)-1;
  Initialize_Params();
  
//   CheckSequence(sequence);
//   S0=encode_seq(sequence+1);
//   seqlen=strlen(sequence)-1;
//   Initialize_Params();
//   make_pair_matrix();//needed for pair matching
//   kT = (temperature+K0)*GASCONST/1000.0;
//   //printf("%.15f\n",kT);
//   ML_base=(double)P->MLbase/100;
//   ML_close=(double)P->MLclosing/100;
//   HPMAX=Allocate2DMatrix( seqlen+1, seqlen+1);
//   for(i=1;i<seqlen+1;++i)
//     {for (j=1;j<seqlen+1;++j)
//         HPMAX[i][j]=0;
//     }
//   MaxNumHairpin(sequence,HPMAX);
//   if(H_MAX>HPMAX[1][seqlen])
// H_MAX=HPMAX[1][seqlen];
//   for (i=0;i<H_MAX+1;++i)
//     {HP[i]=0;
//       HN[i]=0;}
//   HairpinPartition(HP,HN,H_MAX,sequence);
//   double** McCaskillZ;
//   McCaskillZ=runMcCaskill(sequence);
  
  return 0;
}

double** runMcCaskill(char sequence[MAXSIZE]) {
  int i, j, d;
  double **Z;
  double **ZB;
  double **ZM;
  Z  = Allocate2DMatrix(seqlen + 1, seqlen + 1);
  ZB = Allocate2DMatrix(seqlen + 1, seqlen + 1);
  ZM = Allocate2DMatrix(seqlen + 1, seqlen + 1);
  
  for (d = 0; d <= MIN_PAIR_DIST; ++d) {
    for (i = 1; i <= seqlen - d; ++i) {
      SET_Z(i, i + d, 1)
    }
  }
  
  for (d = MIN_PAIR_DIST + 1; d < seqlen; ++d) {
    for (i = 1; i <= seqlen - d; ++i) {
      j = i + d;
      
        if (canBasePair(i, j, sequence)) {
          solveZB(i, j, sequence, ZB, ZM);
        }
        
        solveZM(i, j, sequence, ZB, ZM);
      
      solveZ(i, j, sequence, Z, ZB);
      }
  }
  
  return Z;
}

void solveZ(int i, int j, char sequence[MAXSIZE], double **Z, double **ZB) { 
  int k;
  
  if(j - i < MIN_PAIR_DIST + 1) {
    SET_Z(i, j, 1)
  } else {
    Z[i][j] += Z[i][j - 1];
    Z[j][i] += Z[j - 1][i];
    
    for (k = i; k <= j - MIN_PAIR_DIST - 1; ++k) { 
      // (k, j) is the rightmost base pair in (i, j).
      if (canBasePair(k, j, sequence)) {
        if (k == i) {
          Z[i][j] += ZB[k][j] * exp(-AU_Penalty(i, j, S0) / kT);
        Z[j][i] += ZB[j][k];
        } else {
          Z[i][j] += Z[i][k - 1] * ZB[k][j] * exp(-AU_Penalty(k, j, S0) / kT);
          Z[j][i] += Z[k - 1][i] * ZB[j][k];
        }
      }
    }
  }
}

void solveZB(int i, int j, char sequence[MAXSIZE], double **ZB, double **ZM) { 
  // (i, j) assumed to b.p. in here.
  int k, l;
  
  // In a hairpin, (i + 1, j - 1) all unpaired.
  ZB[i][j] += exp(-HP_Energy(i, j, S0, sequence + 1) / kT);
  ZB[j][i] += 1;
  
  // Interior loop / bulge / stack / multiloop.
  for (k = i + 1; k <= j - MIN_PAIR_DIST - 1; ++k) {
    for (l = max(k + MIN_PAIR_DIST + 1, j - MAX_INTERIOR_DIST - 1); l < j; ++l) {
      if (canBasePair(k, l, sequence)) {
        // In interior loop / bulge / stack with (i, j) and (k, l), (i + 1, k - 1) and (l + 1, j - 1)
        // are all unpaired.
        ZB[i][j] += ZB[k][l] * exp(-IL_Energy(i, j, k, l, S0) / kT);
        ZB[j][i] += ZB[l][k];
          
        // If (i, j) is the closing b.p. of a multiloop, and (k, l) is the rightmost base pair, 
        // there is at least one hairpin between (i + 1, k - 1).
        ZB[i][j] += exp(-(ML_close + MLbasepairAndAUpenalty(j, i, S0)) / kT) * ZB[k][l] * ZM[i + 1][k - 1];
        ZB[j][i] += ZB[l][k] * ZM[k - 1][i + 1];
      }
    }
  }
}

void solveZM(int i, int j, char sequence[MAXSIZE], double **ZB, double **ZM) { 
  int k;
  
  ZM[i][j] += ZM[i][j - 1] * exp(-1 / kT);
  ZM[j][i] += ZM[j - 1][i];
  
  for (k = i; k <= j - MIN_PAIR_DIST - 1; ++k) {
    if (canBasePair(k, j, sequence)) {
      // Only one stem.
      ZM[i][j] += ZB[k][j] * exp(-ML_base * (k - i) / kT);
      ZM[j][i] += ZB[j][k];
      
      // k needs to be greater than MIN_PAIR_DIST + 2 from i to fit more than one stem.
      if (k > i) {
        ZM[i][j] += ZB[k][j] * ZM[i][k - 1] * exp(-ML_base / kT);
        ZM[j][i] += ZB[j][k] * ZM[k - 1][i];
      }
    }
  }
}