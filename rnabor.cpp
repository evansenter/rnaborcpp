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

typedef struct {
  int id;
  int stack[NBPAIRS+1][NBPAIRS+1];
  int hairpin[31];
  int bulge[MAXLOOP+1];
  int internal_loop[MAXLOOP+1];
  int mismatchI[NBPAIRS+1][5][5];
  int mismatchH[NBPAIRS+1][5][5];
  int mismatchM[NBPAIRS+1][5][5];
  int dangle5[NBPAIRS+1][5];
  int dangle3[NBPAIRS+1][5];
  int int11[NBPAIRS+1][NBPAIRS+1][5][5];
  int int21[NBPAIRS+1][NBPAIRS+1][5][5][5];
  int int22[NBPAIRS+1][NBPAIRS+1][5][5][5][5];
  int F_ninio[5];
  double lxc;
  int MLbase;
  int MLintern[NBPAIRS+1];
  int MLclosing;
  int TerminalAU;
  int DuplexInit;
  int TETRA_ENERGY[200];
  char Tetraloops[1401];
  int Triloop_E[40];
  char Triloops[241];
  double temperature;
}  paramT;

paramT *P, p;

int CheckSequence(char sequence[MAXSIZE]);
void *space(unsigned size);
int encode_char(char c);
short *encode_seq(const char *seq);
double **Allocate2DMatrix(int m, int n);
double HP_Energy(int i, int j, short *S0, char* sequence);
double IL_Energy(int i, int j, int ip, int jp, short *S0);
double AU_Penalty(int i, int j, short *S0);
double MLbasepairAndAUpenalty(int i, int j, short *S0);
int canBasePair(int i, int j, char sequence[MAXSIZE]);
double** runMcCaskill(char sequence[MAXSIZE]);
void solveZ(int i, int j, char sequence[MAXSIZE], double **Z, double **ZB);
void solveZB(int i, int j, char sequence[MAXSIZE], double **ZB, double **ZM);
void solveZM(int i, int j, char sequence[MAXSIZE], double **ZB, double **ZM);
void Initialize_Params();
paramT *scale_parameters(void);

paramT *scale_parameters(void)
{
  unsigned int i,j,k,l;
  double tempf;
  /* if ((fabs(p.temperature - temperature)<1e-6)&&(id == p.id)) return &p; */
  tempf = ((temperature+K0)/Tmeasure);
  for (i=0; i<31; i++) 
    p.hairpin[i] = (int) hairpin37[i]*(tempf);
  for (i=0; i<=MIN2(30,MAXLOOP); i++) {
    p.bulge[i] = (int) bulge37[i]*tempf;
    p.internal_loop[i]= (int) internal_loop37[i]*tempf;
  }
  p.lxc = lxc37*tempf;
  for (; i<=MAXLOOP; i++) {
    p.bulge[i] = p.bulge[30]+(int)(p.lxc*log((double)(i)/30.));
    p.internal_loop[i] = p.internal_loop[30]+(int)(p.lxc*log((double)(i)/30.));
  }
  for (i=0; i<5; i++)
    p.F_ninio[i] = (int) F_ninio37[i]*tempf;
   
  for (i=0; (i*7)<strlen(Tetraloops); i++) 
    p.TETRA_ENERGY[i] = TETRA_ENTH37 - (TETRA_ENTH37-TETRA_ENERGY37[i])*tempf;
  for (i=0; (i*5)<strlen(Triloops); i++) 
    p.Triloop_E[i] =  Triloop_E37[i];
   
  p.MLbase = ML_BASE37*tempf;
  for (i=0; i<=NBPAIRS; i++) { /* includes AU penalty */
    p.MLintern[i] = ML_intern37*tempf;
    p.MLintern[i] +=  (i>2)?TerminalAU:0;
  }
  p.MLclosing = ML_closing37*tempf;

  p.TerminalAU = TerminalAU;
  
  p.DuplexInit = DuplexInit*tempf;

  /* stacks    G(T) = H - [H - G(T0)]*T/T0 */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      p.stack[i][j] = enthalpies[i][j] -
	(enthalpies[i][j] - stack37[i][j])*tempf;

  /* mismatches */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<5; j++)
      for (k=0; k<5; k++) {
	p.mismatchI[i][j][k] = mism_H[i][j][k] -
	  (mism_H[i][j][k] - mismatchI37[i][j][k])*tempf;
	p.mismatchH[i][j][k] = mism_H[i][j][k] -
	  (mism_H[i][j][k] - mismatchH37[i][j][k])*tempf;
	p.mismatchM[i][j][k] = mism_H[i][j][k] -
	  (mism_H[i][j][k] - mismatchM37[i][j][k])*tempf;
      }
   
  /* dangles */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<5; j++) {
      int dd;
      dd = dangle5_H[i][j] - (dangle5_H[i][j] - dangle5_37[i][j])*tempf; 
      p.dangle5[i][j] = (dd>0) ? 0 : dd;  /* must be <= 0 */
      dd = dangle3_H[i][j] - (dangle3_H[i][j] - dangle3_37[i][j])*tempf;
      p.dangle3[i][j] = (dd>0) ? 0 : dd;  /* must be <= 0 */
    }
  /* interior 1x1 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
	for (l=0; l<5; l++) 
	  p.int11[i][j][k][l] = int11_H[i][j][k][l] -
	    (int11_H[i][j][k][l] - int11_37[i][j][k][l])*tempf;

  /* interior 2x1 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
	for (l=0; l<5; l++) {
	  int m;
	  for (m=0; m<5; m++)
	    p.int21[i][j][k][l][m] = int21_H[i][j][k][l][m] -
	      (int21_H[i][j][k][l][m] - int21_37[i][j][k][l][m])*tempf;
	}
  /* interior 2x2 loops */
  for (i=0; i<=NBPAIRS; i++)
    for (j=0; j<=NBPAIRS; j++)
      for (k=0; k<5; k++)
	for (l=0; l<5; l++) {
	  int m,n;
	  for (m=0; m<5; m++)
	    for (n=0; n<5; n++)	     
	      p.int22[i][j][k][l][m][n] = int22_H[i][j][k][l][m][n] -
		(int22_H[i][j][k][l][m][n]-int22_37[i][j][k][l][m][n])*tempf;
	}

  strncpy(p.Tetraloops, Tetraloops, 1400);
  strncpy(p.Triloops, Triloops, 240);

  p.temperature = temperature;
  p.id = ++id;
  return &p;
}

inline int LoopEnergy(int n1, int n2, int type, int type_2,
		      int si1, int sj1, int sp1, int sq1) {
  /* compute energy of degree 2 loop (stack bulge or interior) */
  int nl, ns, energy;
  //  printf("in LoopEnergy, n1, n2, type, type_2, si1, sj1, sp1, sq1 are all\n");
  //  printf("%d,%d,%d,%d,%d,%d,%d,%d.\n", n1, n2, type, type_2, si1, sj1,
  // sp1, sq1);
  if (n1>n2) { nl=n1; ns=n2;}
  else {nl=n2; ns=n1;}
  //  printf("nl, ns are %d, %d\n", nl, ns);
  if (nl == 0){
    //  printf("got a stack here!\n");
    return P->stack[type][type_2];    /* stack */
  }
  // printf("no stack?\n");
  if (ns==0) {                       /* bulge */
    energy = (nl<=MAXLOOP)?P->bulge[nl]:
      (P->bulge[30]+(int)(P->lxc*log(nl/30.)));
    if (nl==1) energy += P->stack[type][type_2];
    else {
      if (type>2) energy += P->TerminalAU;
      if (type_2>2) energy += P->TerminalAU;
    }
    return energy;
  }
  else {                             /* interior loop */
    if (ns==1) {
      if (nl==1)                     /* 1x1 loop */
	return P->int11[type][type_2][si1][sj1];
      if (nl==2) {                   /* 2x1 loop */
	if (n1==1)
	  energy = P->int21[type][type_2][si1][sq1][sj1];
	else
	  energy = P->int21[type_2][type][sq1][si1][sp1];
	return energy;
      }
    }
    else if (n1==2 && n2==2)         /* 2x2 loop */
      return P->int22[type][type_2][si1][sp1][sq1][sj1];
    { /* generic interior loop (no else here!)*/
      energy = (n1+n2<=MAXLOOP)?(P->internal_loop[n1+n2]):
	(P->internal_loop[30]+(int)(P->lxc*log((n1+n2)/30.)));

      energy += MIN2(MAX_NINIO, (nl-ns)*P->F_ninio[2]);

      energy += P->mismatchI[type][si1][sj1]+
	P->mismatchI[type_2][sq1][sp1];
    }
  }
  return energy;
}

inline int HairpinE(int size, int type, int si1, int sj1, const char *string) {
  int energy;
  energy = (size <= 30) ? P->hairpin[size] :
    P->hairpin[30]+(int)(P->lxc*log((double) ((size)/30.) ));
  if (tetra_loop)
    if (size == 4) { /* check for tetraloop bonus */
      char tl[7]={0}, *ts;
      strncpy(tl, string, 6);
      if ((ts=strstr(P->Tetraloops, tl)))
	energy += P->TETRA_ENERGY[(ts - P->Tetraloops)/7];
    }
  if (size == 3) {
    char tl[6]={0,0,0,0,0,0}, *ts;
    strncpy(tl, string, 5);
    if ((ts=strstr(P->Triloops, tl)))
      energy += P->Triloop_E[(ts - P->Triloops)/6];

    if (type>2)  /* neither CG nor GC */
      energy += P->TerminalAU; /* penalty for closing AU GU pair */
  }
  else  /* no mismatches for tri-loops */
    energy += P->mismatchH[type][si1][sj1];

  return energy;
}

void *space(unsigned size) {
  void *pointer;
  
  if ( (pointer = (void *) calloc(1, (size_t) size)) == NULL) {
      printf("SPACE allocation failure -> no memory");
  }
  return  pointer;
}

int encode_char(char c) {
  /* return numerical representation of base used e.g. in pair[][] */
  int code;
  if (energy_set>0) code = (int) (c-'A')+1;
  else {
    char *pos;
    pos = strchr(Law_and_Order, c);
    if (pos==NULL) code=0;
    else code = (int) (pos-Law_and_Order);
    if (code>4) code--; /* make T and U equivalent */
  }
  return code;
}

short *encode_seq(const char *seq) {
  unsigned int k,l;
  short *S0_out;
  int i;
  l = strlen(seq);
  S0_out = (short *) space(sizeof(short)*(l+2));
  S0_out[0]=l;
  for (k=1; k<=l; k++) { /* make numerical encoding of seq */
    S0_out[k]= (short) encode_char(toupper(seq[k-1]));
  }
  return S0_out;
}

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

double HP_Energy(int i, int j, short *S0, char* sequence){
  int type, energy_int;

  type=pair[S0[i]][S0[j]];
  energy_int=HairpinE(j-i-1, type, S0[i+1], S0[j-1], sequence+i-1);
  return ((double) energy_int)/100.;
}

double IL_Energy(int i, int j, int ip, int jp, short *S0){
  int type, type2, energy_int, k, l;

  type=pair[S0[i]][S0[j]];
  type2=pair[S0[jp]][S0[ip]];
  energy_int=LoopEnergy(ip-i-1, j-jp-1, type, type2,
			S0[i+1], S0[j-1], S0[ip-1], S0[jp+1]);
  return ((double) energy_int)/100.;
}

double AU_Penalty(int i, int j, short *S0){
  int type, energy_int;
  type=pair[S0[i]][S0[j]];
  energy_int=P->MLintern[type]-P->MLintern[1]; /* 0 or AU penalty */
  return (double) energy_int/100.;
}

double MLbasepairAndAUpenalty(int i, int j, short *S0){
  int type, energy_int;
  int k;
  type=pair[S0[i]][S0[j]];
  energy_int=P->MLintern[type]; 
  return (double) energy_int/100.;
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