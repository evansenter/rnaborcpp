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

double** runMcCaskill(char sequence[MAXSIZE]) {
  int i, j, d;
  double **Z;
  double **ZB;
  double **ZM;
  Z  = Allocate2DMatrix(seqlen + 1, seqlen + 1);
  ZB = Allocate2DMatrix(seqlen + 1, seqlen + 1);
  ZM = Allocate2DMatrix(seqlen + 1, seqlen + 1);
  
  for (i = 0; i < seqlen + 1; ++i) {
    for (j = 0; j < seqlen + 1; ++j) {
  	  Z[i][j]  = 0;
  	  Z[j][i]  = 0;
  	  ZB[i][j] = 0;
  	  ZM[i][j] = 0;
	  }
  }
  
  for (d = MIN_PAIR_DIST + 1; d < seqlen; ++d) {
    for (i = 1; i <= seqlen - d; ++i) {
      j = i + d;
	    
	    if (BP(i, j, sequence)) {
	      solveZB(i, j, sequence, ZB, ZM);
	    }
	    
	    solveZM(i, j, sequence, ZB, ZM);
	  }
  }
  
  for (d = 0; d < seqlen; ++d) {
    for (i = 1; i <= seqlen - d; ++i) {
	    j = i + d;
	    
	    solveZ(i, j, sequence, Z, ZB);
	  }
  }
  
  return Z;
}

int solveZ(int i, int j, char sequence[MAXSIZE], double **Z, double **ZB) { 
  int k;
  
  if(j - i < MIN_PAIR_DIST + 1) {
    Z[i][j] = 1;
    Z[j][i] = 1;
  } else {
    Z[i][j] += Z[i][j - 1];
    Z[j][i] += Z[j - 1][i];
    
    for (k = i; k < j - MIN_PAIR_DIST; ++k) { 
      // (k, j) is the rightmost base pair in (i, j).
	    if (BP(k, j, sequence)) {
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

int solveZB(int i, int j, char sequence[MAXSIZE], double **ZB, double **ZM) { 
  // (i, j) assumed to b.p. in here.
  int k, l;
  
  // In a hairpin, (i + 1, j - 1) all unpaired.
  ZB[i][j] += exp(-HP_Energy(i, j, S0, sequence + 1) / kT);
  ZB[j][i] += 1.0;
  
  // Interior loop / bulge / stack / multiloop.
  for (k = i + 1; k < j - MIN_PAIR_DIST; ++k) {
    for (l = max(k + MIN_PAIR_DIST + 1, j - MAX_INTERIOR_DIST - 1; l < j; ++l) {
      if (BP(k, l, sequence)) {
        // In interior loop / bulge / stack with (i, j) and (k, l), (i + 1, k - 1) and (l + 1, j - 1)
        // are all unpaired.
	      ZB[i][j] += ZB[k][l] * exp(-IL_Energy(i, j, k, l, S0) / kT);
	      ZB[j][i] += ZB[l][k];
	      
	      // If (i, j) is the closing b.p. of a multiloop, and (k, l) is the rightmost base pair, 
        // there is at least one hairpin between (i + 1, k - 1) in order to be a multiloop, so
        // k needs to be more than (MIN_PAIR_DIST + 2) from i to have room for the branch.
  	    if (k > i + MIN_PAIR_DIST + 2) {
          ZB[i][j] += exp(-(ML_close + MLbasepairAndAUpenalty(j, i, S0)) / kT) * ZB[k][l] * ZM[i + 1][k - 1];
          ZB[j][i] += ZB[l][k] * ZM[k - 1][i + 1];
  	    }
	    }
	  }
  }
}

int solveZM(int i, int j, char sequence[MAXSIZE], double **ZB, double **ZM) { 
  int k;
  
  for (k = i; k < j - MIN_PAIR_DIST; ++k) {
    if (BP(k, j, sequence)) {
      // Only one stem.
      ZM[i][j] += ZB[k][j] * exp(-ML_base * (k - i) / kT);
      ZM[j][i] += ZB[j][k];
      
      // k needs to be greater than MIN_PAIR_DIST + 2 from i to fit more than one stem.
      if (k > i + MIN_PAIR_DIST + 2) {
        ZM[i][j] += ZB[k][j] * ZM[i][k - 1] * exp(-ML_base / kT);
        ZM[j][i] += ZB[j][k] * ZM[k - 1][i];
      }
    }
  }
}
