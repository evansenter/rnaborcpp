#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "myConst.h"
#include "fold.h"
#include "energy_const.h"
#include "fold_vars.h"
#include "pair_mat.h"
#include "convert_Vienna.h"
#include "params.h"
#include <limits.h>
#include "RNAhairpin.h"
#include "misc.h"
#define MIN_PAIR_DIST 3
#define MAX_INTERIOR_DIST 30

double** runMcCaskill(char *sequence) {
  int x, i, j, d;
  double **complexRoots;
  double ***Z;
  double ***ZB;
  double ***ZM;
  char structure[seqlen + 1];
  
  complexRoots = Allocate2DMatrix(seqlen + 1, 4);
  Z            = Allocate3DMatrix(seqlen + 1, seqlen + 1, 2);
  ZB           = Allocate3DMatrix(seqlen + 1, seqlen + 1, 2);
  ZM           = Allocate3DMatrix(seqlen + 1, seqlen + 1, 2);
  
  // For now the default structure is the unpaired structure, this can easily be changed later on.
  for (i = 0; i <= seqlen; ++i) {
    structure[i] = '.';
  }
  
  for (x = 0; x < seqlen + 1; ++x) {
    // Set this iteration's root of unity values.
    complexRoots[x][0] = cos(2 * M_PI * x / (seqlen + 1));
    complexRoots[x][1] = sin(2 * M_PI * x / (seqlen + 1));
    complexRoots[x][2] = 0;
    
    // Reinitialize the 3 matricies.
    for (i = 0; i < seqlen + 1; ++i) {
      for (j = 0; j < seqlen + 1; ++j) {
    	  Z[i][j]  = 0;
    	  ZB[i][j] = 0;
    	  ZM[i][j] = 0;
  	  }
    }

    // Populate ZB / ZM
    for (d = MIN_PAIR_DIST + 1; d < seqlen; ++d) {
      for (i = 1; i <= seqlen - d; ++i) {
        j = i + d;

  	    if (BP(i, j, sequence)) {
  	      solveZB(i, j, x, sequence, structure, ZB, ZM, complexRoots);
  	    }

  	    solveZM(i, j, x, sequence, structure, ZB, ZM, complexRoots);
  	  }
    }

    // Populate Z
    for (d = 0; d < seqlen; ++d) {
      for (i = 1; i <= seqlen - d; ++i) {
  	    j = i + d;

  	    solveZ(i, j, x, sequence, structure, Z, ZB, complexRoots);
  	  }
    }
    
    // Save the complex solution for McCaskill's recursions evaluated at (complexRoots[x][0] + complexRoots[x][1] * i)
    complexRoots[x][2] = Z[1][seqlen][0];
    complexRoots[x][3] = Z[1][seqlen][1];
  }
  
  return complexRoots;
}

int solveZ(int i, int j, int x, char *sequence, char *structure, double ***Z, double ***ZB, double **complexRoots) { 
  int k;
  
  if(j - i < MIN_PAIR_DIST + 1) {
    Z[i][j] = 1;
    Z[j][i] = 1;
  } else {
    Z[i][j] += Z[i][j - 1] * pow(xi, pairedIn(i, j, j, structure) ? 1 : 0);
    Z[j][i] += Z[j - 1][i];
    
    for (k = i; k <= j - MIN_PAIR_DIST - 1; ++k) { 
      // (k, j) is the rightmost base pair in (i, j).
	    if (BP(k, j, sequence)) {
	      if (k == i) {
		      Z[i][j] += ZB[k][j] * exp(-AU_Penalty(i, j, S0) / kT) * pow(xi, bpDifference(i, j, structure, 0, 0, 0, 0, k, j));
		      Z[j][i] += ZB[j][k];
		    } else {
		      Z[i][j] += Z[i][k - 1] * ZB[k][j] * exp(-AU_Penalty(k, j, S0) / kT) * pow(xi, bpDifference(i, j, structure, 0, 0, i, k - 1, k, j));
		      Z[j][i] += Z[k - 1][i] * ZB[j][k];
		    }
	    }
	  }
  }
}

int solveZB(int i, int j, int x, char *sequence, char *structure, double ***ZB, double ***ZM, double **complexRoots) { 
  // (i, j) assumed to b.p. in here.
  int k, l;
  
  // In a hairpin, (i + 1, j - 1) all unpaired.
  ZB[i][j] += exp(-HP_Energy(i, j, S0, sequence + 1) / kT) * pow(xi, bpDifference(i, j, structure, i, j, 0, 0, 0, 0));
  ZB[j][i] += 1;
  
  // Interior loop / bulge / stack / multiloop.
  for (k = i + 1; k <= j - MIN_PAIR_DIST - 1; ++k) {
    for (l = max(k + MIN_PAIR_DIST + 1, j - MAX_INTERIOR_DIST - 1); l < j; ++l) {
      if (BP(k, l, sequence)) {
        // In interior loop / bulge / stack with (i, j) and (k, l), (i + 1, k - 1) and (l + 1, j - 1)
        // are all unpaired.
	      ZB[i][j] += ZB[k][l] * exp(-IL_Energy(i, j, k, l, S0) / kT) * pow(xi, bpDifference(i, j, structure, i, j, k, l, 0, 0));
	      ZB[j][i] += ZB[l][k];
	      
	      // If (i, j) is the closing b.p. of a multiloop, and (k, l) is the rightmost base pair, 
        // there is at least one hairpin between (i + 1, k - 1).
        ZB[i][j] += exp(-(ML_close + MLbasepairAndAUpenalty(j, i, S0)) / kT) * ZB[k][l] * ZM[i + 1][k - 1] * pow(xi, bpDifference(i, j, structure, i, j, i + 1, k - 1, k, l));
        ZB[j][i] += ZB[l][k] * ZM[k - 1][i + 1];
	    }
	  }
  }
}

int solveZM(int i, int j, int x, char *sequence, char *structure, double ***ZB, double ***ZM, double **complexRoots) { 
  int k;
  
  ZM[i][j] += ZM[i][j - 1] * exp(-1 / kT) * pow(xi, pairedIn(i, j, j, structure) ? 1 : 0);
  ZM[j][i] += ZM[j - 1][i];
  
  for (k = i; k <= j - MIN_PAIR_DIST - 1; ++k) {
    if (BP(k, j, sequence)) {
      // Only one stem.
      ZM[i][j] += ZB[k][j] * exp(-ML_base * (k - i) / kT) * pow(xi, bpDifference(i, j, structure, 0, 0, k, j, 0, 0));
      ZM[j][i] += ZB[j][k];
      
      // k needs to be greater than MIN_PAIR_DIST + 2 from i to fit more than one stem.
      if (k > i) {
        ZM[i][j] += ZB[k][j] * ZM[i][k - 1] * exp(-ML_base / kT) * pow(xi, bpDifference(i, j, structure, 0, 0, i, k - 1, k, j));
        ZM[j][i] += ZB[j][k] * ZM[k - 1][i];
      }
    }
  }
}

double* addComplex(double *complexA, double *complexB, double *complexSum) {
  complexSum[0] = complexA[0] + complexB[0];
  complexSum[1] = complexA[1] + complexB[1];
  return complexSum;
}

double* scalarComplex(double *complex, double scalar, double *complexScalar) {
  complexScalar[0] = complex[0] * scalar;
  complexScalar[1] = complex[1] * scalar;
  return complexScalar;
}

double* productComplex(double *complexA, double *complexB, double *complexProduct) {
  complexProduct[0] = complexA[0] * complexB[0] - complexA[1] * complexB[1];
  complexProduct[1] = complexA[0] * complexB[1] + complexA[1] * complexB[0];
  return complexProduct;
}

double* powerComplex(double *complex, int power, double *complexPower) {
  // This only works for integer powers of roots of unity, because we assume the radius of the complex number in polar form
  // is 1, and pow(1, x) == 1
  return rectangularRootComplex(arctan(complex[1] / complex[0]) * power, complexPower);
}

double* rectangularRootComplex(double angle, double *complexResult) {
  complexResult[0] = cos(angle);
  complexResult[1] = sin(angle);
  return complexResult;
}