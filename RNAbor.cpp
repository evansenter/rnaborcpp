#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h> 
#include "fold.h"
#include "energy_const.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "convert_Vienna.h"
#include "params.h"
#include "myConst.h"
#include "misc.h"
#include "RNAbor.h"
#include "McCaskill.h"
#include "limits.h"

double kT;
double ML_base;
double ML_close;
int    seqlen;
short  *S0;
double SCALING_FACTOR;
int    PRECISION;

PRIVATE short *encode_seq(const char *seq);

int main(int argc, char *argv[]){
  char sequence[MAXSIZE], structure[MAXSIZE];
  int i, j;
  
  if (argc == 9) {
    i = 1;
    while (i < argc) {
  	  if (!strcmp(argv[i], "-s")) {
        sequence[0] = '@';
        strncpy(sequence + 1, argv[i + 1], strlen(argv[i + 1])); // The sequence starts from index 1
        sequence[strlen(argv[i + 1]) + 1] = '\0';                // Remember to end the sequence
        i += 2;
      } else if (!strcmp(argv[i], "-r")) {
        strncpy(structure, argv[i + 1], strlen(argv[i + 1])); // The structure starts from index 0
        structure[strlen(argv[i + 1])] = '\0';                // Remember to end the structure
        i += 2;
      } else if (!strcmp(argv[i], "-c")) {
        SCALING_FACTOR = atof(argv[i + 1]);
        i += 2;
      } else if (!strcmp(argv[i], "-p")) {
        PRECISION = atoi(argv[i + 1]);
        i += 2;
      } else {
        printf("Error: Unrecognizable Flag!\n");
        printf("Usage:");
        printf("%s -s SEQUENCE -r STRUCTURE -c SCALING -p PRECISION\n", argv[0]);
        printf("SEQUENCE is the RNA sequence.\n");
        printf("STRUCTURE is the RNA structure.\n");
        printf("SCALING is used internally as the scaling factor for the recursions.\n");
        printf("PRECISION dictates the number of signigicant digits in the resulting distribution.\n");
        printf("The length of RNA sequence should be less than %d and larger than or equal to 5.\n", MAXSIZE);
        exit(1);
      }
  	}
    CheckSequenceAndStructure(sequence, structure);
    S0     = encode_seq(sequence + 1);
    seqlen = strlen(sequence) - 1;
    
    Initialize_Params();
    make_pair_matrix(); // Needed for pair matching.
    
    kT       = (temperature + K0) * GASCONST / 1000.0;
    ML_base  = (double)P -> MLbase / 100;
    ML_close = (double)P -> MLclosing / 100;
    std::complex<double>** McCaskillZ;
    
    McCaskillZ = runMcCaskill(sequence, structure);
  } else {
    printf("Usage:");
    printf("%s -s SEQUENCE -r STRUCTURE -c SCALING -p PRECISION\n", argv[0]);
    printf("SEQUENCE is the RNA sequence.\n");
    printf("STRUCTURE is the RNA structure.\n");
    printf("SCALING is used internally as the scaling factor for the recursions.\n");
    printf("PRECISION dictates the number of signigicant digits in the resulting distribution.\n");
    printf("The length of RNA sequence should be less than %d and larger than or equal to 5.\n", MAXSIZE);
    exit(1);
  }
  return 0;
}

/* Below modified from fold.c */
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
