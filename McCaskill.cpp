// GGGGGCCCCCGGGGGCCCCCGGGGGCCCCC
// 19049760

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
#include "McCaskill.h"
#include <lapackpp.h>
#define MIN_PAIR_DIST 3
#define MAX_INTERIOR_DIST 30
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define ZERO_C dcomplex(0, 0)
#define ONE_C  dcomplex(1, 0)
#define SET_Z(i, j, value) \
Z[i][j] = value; \
Z[j][i] = value;

dcomplex** runMcCaskill(char sequence[MAXSIZE]) {
  // Variable declarations.
  int root, i, j, d, *basePairs, **basePairCounts;
  char structure[seqlen];
  
  dcomplex **Z           = new dcomplex*[seqlen + 1];
  dcomplex **ZB          = new dcomplex*[seqlen + 1];
  dcomplex **ZM          = new dcomplex*[seqlen + 1];
  dcomplex **rootsOfUnity = new dcomplex*[seqlen + 1];
  
  // Matrix allocation.
  for (i = 0; i <= seqlen; ++i) {
    Z[i]              = new dcomplex[seqlen + 1];
    ZB[i]             = new dcomplex[seqlen + 1];
    ZM[i]             = new dcomplex[seqlen + 1];
    rootsOfUnity[i]    = new dcomplex[2];
    rootsOfUnity[i][0] = dcomplex(cos(2 * M_PI * i / (seqlen + 1)), sin(2 * M_PI * i / (seqlen + 1)));
  }
  
  // This will need to be parameterized.
  for (i = 0; i < seqlen; ++i) {
    structure[i] = '.';
    structure[i] = seqlen % 2 != 0 && i + 1 == seqlen ? '.' : (i % 2 == 0 ? '(' : ')');
  }
  
  basePairs      = getBasePairList(structure);
  basePairCounts = fillBasePairCounts(basePairs, seqlen);
  
  // Start main recursions.
  for (root = 0; root <= seqlen; ++root) {
    // Flush the matrices.
    for (i = 0; i <= seqlen; ++i) {
      for (j = 0; j <= seqlen; ++j) {
        Z[i][j]  = ZERO_C;
        ZB[i][j] = ZERO_C;
        ZM[i][j] = ZERO_C;
      }
    }
  
    // Set base case for Z.
    for (d = 0; d <= MIN_PAIR_DIST; ++d) {
      for (i = 1; i <= seqlen - d; ++i) {
        SET_Z(i, i + d, ONE_C)
      }
    }
    
    for (d = MIN_PAIR_DIST + 1; d < seqlen; ++d) {
      for (i = 1; i <= seqlen - d; ++i) {
        j = i + d;
      
        if (BP(i, j, sequence)) {
          solveZB(i, j, rootsOfUnity[root][0], sequence, basePairs, basePairCounts, ZB, ZM);
        }
        
        solveZM(i, j, rootsOfUnity[root][0], sequence, basePairs, basePairCounts, ZB, ZM);
      
        solveZ(i, j, rootsOfUnity[root][0], sequence, basePairs, basePairCounts, Z, ZB);
      }
    }
    
    rootsOfUnity[root][1] = Z[1][seqlen];
  }
  
  // LaGenMatComplex A(row, col);
  // LaVectorComplex X(col);
  // LaVectorComplex B(col);
  // 
  // for (i = 0; i < row; i++) {
  //   for (j = 0; j < col; j++) {
  //     A(i, j).r = i + j + 0.5;
  //     A(i, j).i = -(i + j + 0.5);
  //   }
  //   
  //   B(i).r = i + 0.5;
  //   B(i).i = -(i + 0.5);
  // }
  // 
  // LaLinearSolveIP(A, X, B);
  
  return Z;
}

void solveZ(int i, int j, dcomplex x, char sequence[MAXSIZE], int *basePairs, int **basePairCounts, std::complex<double> **Z, std::complex<double> **ZB) { 
  int k;
  
  if(j - i < MIN_PAIR_DIST + 1) {
    SET_Z(i, j, 1)
  } else {
    Z[i][j] += Z[i][j - 1] * pow(x, basePairs[j] < j && basePairs[j] >= i ? 1 : 0);
    Z[j][i] += Z[j - 1][i];
    
    for (k = i; k <= j - MIN_PAIR_DIST - 1; ++k) { 
      // (k, j) is the rightmost base pair in (i, j).
      if (BP(k, j, sequence)) {
        if (k == i) {
          Z[i][j] += ZB[k][j] * exp(-AU_Penalty(i, j, S0) / kT) * x;
          Z[j][i] += ZB[j][k];
        } else {
          Z[i][j] += Z[i][k - 1] * ZB[k][j] * exp(-AU_Penalty(k, j, S0) / kT) * pow(x, basePairCounts[i][j] - basePairCounts[i][k - 1] - basePairCounts[k][j]);
          Z[j][i] += Z[k - 1][i] * ZB[j][k];
        }
      }
    }
  }
}

void solveZB(int i, int j, dcomplex x, char sequence[MAXSIZE], int *basePairs, int **basePairCounts, std::complex<double> **ZB, std::complex<double> **ZM) { 
  // (i, j) assumed to b.p. in here.
  int k, l;
  
  // In a hairpin, (i + 1, j - 1) all unpaired.
  ZB[i][j] += exp(-HP_Energy(i, j, S0, sequence + 1) / kT) * x;
  ZB[j][i] += 1;
  
  // Interior loop / bulge / stack / multiloop.
  for (k = i + 1; k <= j - MIN_PAIR_DIST - 1; ++k) {
    for (l = MAX(k + MIN_PAIR_DIST + 1, j - MAX_INTERIOR_DIST - 1); l < j; ++l) {
      if (BP(k, l, sequence)) {
        // In interior loop / bulge / stack with (i, j) and (k, l), (i + 1, k - 1) and (l + 1, j - 1)
        // are all unpaired.
        ZB[i][j] += ZB[k][l] * exp(-IL_Energy(i, j, k, l, S0) / kT) * x;
        ZB[j][i] += ZB[l][k];
          
        // If (i, j) is the closing b.p. of a multiloop, and (k, l) is the rightmost base pair, 
        // there is at least one hairpin between (i + 1, k - 1).
        ZB[i][j] += exp(-(ML_close + MLbasepairAndAUpenalty(j, i, S0)) / kT) * ZB[k][l] * ZM[i + 1][k - 1] * x;
        ZB[j][i] += ZB[l][k] * ZM[k - 1][i + 1];
      }
    }
  }
}

void solveZM(int i, int j, dcomplex x, char sequence[MAXSIZE], int *basePairs, int **basePairCounts, std::complex<double> **ZB, std::complex<double> **ZM) { 
  int k;
  
  ZM[i][j] += ZM[i][j - 1] * exp(-1 / kT) * x;
  ZM[j][i] += ZM[j - 1][i];
  
  for (k = i; k <= j - MIN_PAIR_DIST - 1; ++k) {
    if (BP(k, j, sequence)) {
      // Only one stem.
      ZM[i][j] += ZB[k][j] * exp(-ML_base * (k - i) / kT) * x;
      ZM[j][i] += ZB[j][k];
      
      // k needs to be greater than MIN_PAIR_DIST + 2 from i to fit more than one stem.
      if (k > i) {
        ZM[i][j] += ZB[k][j] * ZM[i][k - 1] * exp(-ML_base / kT) * x;
        ZM[j][i] += ZB[j][k] * ZM[k - 1][i];
      }
    }
  }
}

int* getBasePairList(char *structure) {
  // Assumes that structure is 0-indexed.
  int length     = strlen(structure);
  int *stack     = (int *) calloc(length, sizeof(int));
  int *pairsList = (int *) calloc(length + 1, sizeof(int));
  int i, j = 0;
  char symbol;

  // First position holds the number of base pairs.
  for (i = 1; i <= length; ++i) {
    pairsList[i] = -1;
  }

  for (i = 1; i <= length; ++i) {
    symbol = structure[i - 1];
    
    if (symbol == '(') {
      stack[j++] = i;
    } else if (symbol == ')') {
      if (j == 0) {
        std::cout << "There is something wrong with the structure, too many ')'" << std::endl;
        pairsList[0] = -1; 
        return pairsList;
      } else {
        pairsList[stack[--j]] = i;
	      pairsList[i] = stack[j];
	      pairsList[0]++;
      }
    }
  }

  if (j != 0) {
    std::cout << "There is something wrong with the structure, too many '('" << std::endl;
    pairsList[0] = -2;
  }
  
  free(stack);

  return pairsList;
}

/* Number of base pairs in the region i to j in basePairs */
int numberOfBasePairs(int i, int j, int *basePairs) {
  int n = 0;
  int k;
  for (k = i; k <= j; ++k) {
    // If position k opens a b.p. '(' and is closed within the range [i, j], increment n
    if (k < basePairs[k] && basePairs[k] <= j) {
      n++;
    }
  }
  
  return n;
}

int** fillBasePairCounts(int *basePairs, int n) {
  int i, d, **basePairCounts;
  
  basePairCounts = (int **) calloc(n + 1, sizeof(int *));
  for (i = 1; i <= n; ++i) {
    basePairCounts[i] = (int *) calloc(n + 1, sizeof(int));
  }
  
  for (d = MIN_PAIR_DIST + 1; d < n; ++d) {
    for (i = 1; i <= n - d; ++i) {
      basePairCounts[i][i + d] = numberOfBasePairs(i, i + d, basePairs);
    }
  }
  
  return basePairCounts;
}