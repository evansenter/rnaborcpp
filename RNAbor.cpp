#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h> 
#include"fold.h"
#include"energy_const.h"
#include"fold_vars.h"
#include"pair_mat.h"
#include"convert_Vienna.h"
#include"params.h"
#include"myConst.h"
#include"misc.h"
#include"RNAbor.h"
#include"McCaskill.h"
#include"limits.h"

double kT;
double ML_base;
double ML_close;
int seqlen;
short *S0;

PRIVATE short *encode_seq(const char *seq);

int main(int argc, char *argv[]){
  char sequence[MAXSIZE],filename[FILENAMELENGTH];
  FILE *inputfile;
  int i,j,h;
  double **HPMAX;
  double HP[MAXSIZE],HN[MAXSIZE];
  int H_MAX=INT_MAX;
  std::complex<double> totalpar; //total partition function
  if (argc==3||argc==5)
    {
      i=1;
      while(i<argc)
	{
	  if(!strcmp(argv[i], "-s"))
	    {
	      sequence[0]='@';
	      strncpy(sequence+1,argv[i+1],strlen(argv[i+1])); // the sequence starts from index 1
	      sequence[strlen(argv[i+1])+1]='\0'; //remember to end the sequence
	      i+=2;
	    }
	  else if(!strcmp(argv[i], "-n"))
	    {
	      H_MAX=atoi(argv[i+1]);
	      i+=2;
	    }
	  else if(argc==3)
	    {
	     i+=2;
	    }
	  else
	    {
	      printf("Error: Unrecognizable Flag!\n");
	      printf("Usage:");
	      printf("%s -s SEQUENCE [-n NUMBER]\n",argv[0]);
	      printf("SEQUENCE is the RNA sequence.\n");
	      printf("NUMBER is the maximum number of hairpin that one considers.\n");
	      printf("If not specified, the program will calculate the maximum possible number of hairpins given the current RNA sequence, and use that number as an upper bound.\n");
	      printf("The length of RNA sequence should be less than %d and larger than or equal to 5.\n",MAXSIZE);
	      exit(1);
	    }
	}
    CheckSequence(sequence);
    S0=encode_seq(sequence+1);
    seqlen=strlen(sequence)-1;
    Initialize_Params();
    make_pair_matrix();//needed for pair matching
    kT = (temperature+K0)*GASCONST/1000.0;
    //printf("%.15f\n",kT);
    ML_base=0;//(double)P->MLbase/100;
    ML_close=0;//(double)P->MLclosing/100;
    std::complex<double>** McCaskillZ;
    McCaskillZ=runMcCaskill(sequence);
    totalpar=McCaskillZ[1][seqlen];
    printf("The RNA sequence is %s:\n",sequence+1);
    printf("The total number of structures is %.0f.\n", McCaskillZ[seqlen][1].real());
  }
  else{
    printf("Usage:");
    printf("%s -s SEQUENCE [-n NUMBER]\n",argv[0]);
    printf("SEQUENCE is the RNA sequence\n");
    printf("NUMBER is the maximum number of hairpin that one considers.\n");
    printf("If not specified, the program will calculate the maximum possible number of hairpins given the current RNA sequence, and use that number as an upper bound\n");
    printf("The length of RNA sequence should be less than %d and larger than or equal to 5\n",MAXSIZE);
    exit(1);
  }
  return 0;
}

/*Below modified from fold.c*/
PRIVATE short *encode_seq(const char *seq) {
  unsigned int k,l;
  short *S0_out;
  int i;
  l = strlen(seq);
  S0_out = (short *) calloc(1, (size_t) sizeof(short)*(l+2));
  S0_out[0]=l;
  for (k=1; k<=l; k++) { /* make numerical encoding of seq */
    S0_out[k]= (short) encode_char(toupper(seq[k-1]));
  }
  return S0_out;
}
