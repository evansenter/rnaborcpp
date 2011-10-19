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

double** runMcCaskill(char sequence[MAXSIZE]) {
  int i, j, d;
  double **McZ;
  double **McZB;
  double **McZM;
  double **McZM1;
  McZ   = Allocate2DMatrix(seqlen + 1, seqlen + 1);
  McZB  = Allocate2DMatrix(seqlen + 1, seqlen + 1);
  McZM  = Allocate2DMatrix(seqlen + 1, seqlen + 1);
  McZM1 = Allocate2DMatrix(seqlen + 1, seqlen + 1);
  
  for (i = 0; i < seqlen + 1; ++i) {
    for (j = 0; j < seqlen + 1; ++j) {
  	  McZB[i][j]  = 0;
  	  McZ[i][j]   = 0;
  	  McZM[i][j]  = 0;
  	  McZM1[i][j] = 0;
	  }
  }
  
  for (d = 4; d < seqlen; ++d) {
    for (i = 1; i <= seqlen - d; ++i) {
      j = i + d;
	    
	    if (BP(i, j, sequence)) {
	      McGetZB(i, j, sequence, McZM1, McZM, McZB);
	    }
	    
	    McGetZM1(i, j, sequence, McZM1, McZB);
	    McGetZM(i, j, sequence, McZM1, McZM);
	  }
  }
  
  for (d = 0; d < seqlen; ++d) {
    for (i = 1; i <= seqlen - d; ++i) {
	    j = i + d;
	    
	    McGetZ(i, j, sequence, McZ, McZB);
	  }
  }
  
  return McZ;
}

int McGetZ(int i, int j, char sequence[MAXSIZE], double **Z, double **ZB) { 
  int r;
  
  if(j - i < 4) {
    Z[i][j] = exp(0);
    Z[j][i] = 1;
  } else {
    Z[i][j] += Z[i][j - 1];
    Z[j][i] += Z[j - 1][i];
    
    for (r = i; r < j - 3; ++r) { 
	    if (BP(r, j, sequence)) {
	      if (r == i) {
		      Z[i][j] += ZB[r][j] * exp(-AU_Penalty(i, j, S0) / kT);
		      Z[j][i] += ZB[j][r];
		    } else {
		      Z[i][j] += Z[i][r - 1] * ZB[r][j] * exp(-AU_Penalty(r, j, S0) / kT);
		      Z[j][i] += Z[r - 1][i] * ZB[j][r];
		    }
	    }
	  }
  }
}

int McGetZB(int i, int j, char sequence[MAXSIZE], double **ZM1, double **ZM, double **ZB) { 
  int l, r;
  
  ZB[i][j] += exp(-HP_Energy(i, j, S0, sequence + 1) / kT);
  ZB[j][i] += 1.0;
  
  for (l = i + 1; l < min(i + 30, j - 5) + 1; ++l) {
    for (r = max(l + 4, j - (30 - (l - i))); r < j; ++r) {
      if (BP(l, r, sequence)) {
	      ZB[i][j] += ZB[l][r] * exp(-IL_Energy(i, j, l, r, S0) / kT);
	      ZB[j][i] += ZB[r][l];
	    }
	  }
  }
  
  for (r = i + 6; r < j - 4; ++r) {
    ZB[i][j] += exp(-(ML_close + MLbasepairAndAUpenalty(j, i, S0)) / kT) * ZM[i + 1][r - 1] * ZM1[r][j - 1];
    ZB[j][i] += ZM[r - 1][i + 1] * ZM1[j - 1][r];
  }
}

int McGetZM(int i, int j, char sequence[MAXSIZE], double **ZM1, double **ZM) { 
  int r;
  
  for (r = i; r < j - 3; ++r) {
    ZM[i][j] += ZM1[r][j] * exp(-ML_base * (r - i) / kT);
    ZM[j][i] += ZM1[j][r];
  }
     
  for (r = i + 5; r < j - 3; ++r) {
    ZM[i][j] += ZM[i][r-1] * ZM1[r][j];
    ZM[j][i] += ZM[r-1][i] * ZM1[j][r];
  }
}

int McGetZM1(int i, int j, char sequence[MAXSIZE], double **ZM1, double **ZB) { 
  int r;
  
  for (r = i + 4; r < j + 1; ++r) {
    ZM1[i][j] += ZB[i][r] * exp(-(ML_base * (j - r) + MLbasepairAndAUpenalty(i, r, S0)) / kT);
    ZM1[j][i] += ZB[r][i];
  }
}