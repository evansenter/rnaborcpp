// GGGGGCCCCCGGGGGCCCCCGGGGGCCCCC
// 19049760
// 
// CGUUGUGACCGGAAUGGAGUGGUGUCUGGUCUGGGACAUGUAUACCGCAAUACAGCUCCUCUCUUGGUCCUUAGUCCCUGGUUUUGUCACGCUUAAUCCC
// 596353231727330381639733
// 
// GGGGGCCCCCGGGGGCCCCCGGGGGCCCCCGGGGGCCCCCGGGGGCCCCCGGGGGCCCCC
// 11843733310154873 (according to webserver)
// 11540343231278612 (according to Zuker [with MAX_INTERIOR_DIST 30])
// 11540343268028146 (according to Zuker [without MAX_INTERIOR_DIST])
// 11843733310160115 (according to Nussinov)

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
#define STRUCTURE_COUNT 1
#define MIN_PAIR_DIST 3
#define MAX_INTERIOR_DIST 30
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define ZERO_C dcomplex(0.0, 0.0)
#define SET_Z(i, j, value) \
Z[i][j] = value; \
Z[j][i] = value;

dcomplex** runMcCaskill(char sequence[MAXSIZE]) {
  // Variable declarations.
  int root, i, j, k, d, *basePairs, **basePairCounts;
  char structure[seqlen];
  
  dcomplex **Z            = new dcomplex*[seqlen + 1];
  dcomplex **ZB           = new dcomplex*[seqlen + 1];
  dcomplex **ZM           = new dcomplex*[seqlen + 1];
  dcomplex **rootsOfUnity = new dcomplex*[seqlen + 1];
  
  // Roots of unity hack.
  // Matrix allocation.
  for (i = 0; i <= seqlen; ++i) {
    Z[i]               = new dcomplex[seqlen + 1];
    ZB[i]              = new dcomplex[seqlen + 1];
    ZM[i]              = new dcomplex[seqlen + 1];
    rootsOfUnity[i]    = new dcomplex[2];
    rootsOfUnity[i][0] = dcomplex(cos(2 * M_PI * i / (seqlen + 1)), sin(2 * M_PI * i / (seqlen + 1)));
  }
  
  // This will need to be parameterized.
  for (i = 0; i < seqlen; ++i) {
    structure[i] = '.';
  }
  
  basePairs      = getBasePairList(structure);
  basePairCounts = fillBasePairCounts(basePairs, seqlen);
  
  // Start main recursions (root < seqlen / 2 + 1 is an optimization for roots of unity).
  for (root = 0; root < seqlen / 2 + 1; ++root) {
    std::cout << '.' << std::flush;
    
    // Flush the matrices.
    for (i = 0; i <= seqlen; ++i) {
      for (j = 0; j <= seqlen; ++j) {
        Z[i][j]  = ZERO_C;
        ZB[i][j] = ZERO_C;
        ZM[i][j] = ZERO_C;
      }
    }
    
    for (d = 0; d <= MIN_PAIR_DIST; ++d) {
      for (i = 1; i <= seqlen - d; ++i) {
        SET_Z(i, i + d, dcomplex(1, 0))
        // SET_Z(i, i + d, dcomplex(1.0 / (d + 1), 0))
      }
    }
    
    for (d = MIN_PAIR_DIST + 1; d < seqlen; ++d) {
      for (i = 1; i <= seqlen - d; ++i) {
        j = i + d;
      
        // if (BP(i, j, sequence)) {
        //   solveZB(i, j, rootsOfUnity[root][0], sequence, basePairs, basePairCounts, ZB, ZM);
        // }
        // 
        // solveZM(i, j, rootsOfUnity[root][0], sequence, basePairs, basePairCounts, ZB, ZM);
        // 
        // solveZ(i, j, rootsOfUnity[root][0], sequence, basePairs, basePairCounts, Z, ZB);
        
        Z[i][j] = Z[i][j - 1] * pow(rootsOfUnity[root][0], jPairedIn(i, j, basePairs));// / dcomplex(2.5, 0);
        
        for (k = i; k <= j - MIN_PAIR_DIST - 1; ++k) { 
          if (BP(k, j, sequence)) {
            if (k == i) {
              Z[i][j] += Z[k + 1][j - 1] * pow(rootsOfUnity[root][0], basePairCounts[i][j] - basePairCounts[k + 1][j - 1] + jPairedTo(k, j, basePairs));// / dcomplex(6.25, 0);
            } else {
              Z[i][j] += Z[i][k - 1] * Z[k + 1][j - 1] * pow(rootsOfUnity[root][0], basePairCounts[i][j] - basePairCounts[i][k - 1] - basePairCounts[k + 1][j - 1] + jPairedTo(k, j, basePairs));// / dcomplex(6.25, 0);
            }
          }
        }
      }
    }
    
    rootsOfUnity[root][1] = Z[1][seqlen];
  }
  
  // Optimization leveraging complementarity of roots of unity.
  for (i = root - 1; root <= seqlen && i > 0; ++root, --i) {
    rootsOfUnity[root][1] = dcomplex(rootsOfUnity[i][1].real(), -rootsOfUnity[i][1].imag());
  }
  
  std::cout << std::endl;
  
  solveLinearSystem(rootsOfUnity);
  
  return Z;
}

void solveZ(int i, int j, dcomplex x, char sequence[MAXSIZE], int *basePairs, int **basePairCounts, dcomplex **Z, dcomplex **ZB) { 
  int k;
  
  Z[i][j] += Z[i][j - 1] * pow(x, jPairedIn(i, j, basePairs));
    
  #ifdef STRUCTURE_COUNT
    Z[j][i] += Z[j - 1][i];
  #endif
    
  for (k = i; k <= j - MIN_PAIR_DIST - 1; ++k) { 
    // (k, j) is the rightmost base pair in (i, j).
    if (BP(k, j, sequence)) {
      if (k == i) {
        Z[i][j] += ZB[k][j] * pow(x, jPairedIn(i, j, basePairs));
        // Z[i][j] += ZB[k][j] * exp(-AU_Penalty(i, j, S0) / kT) * pow(x, jPairedIn(i, j, basePairs));
          
        #ifdef STRUCTURE_COUNT
          Z[j][i] += ZB[j][k];
        #endif
      } else {
        Z[i][j] += Z[i][k - 1] * ZB[k][j] * pow(x, basePairCounts[i][j] - basePairCounts[i][k - 1] - basePairCounts[k][j]);
        // Z[i][j] += Z[i][k - 1] * ZB[k][j] * exp(-AU_Penalty(k, j, S0) / kT) * pow(x, basePairCounts[i][j] - basePairCounts[i][k - 1] - basePairCounts[k][j]);
          
        #ifdef STRUCTURE_COUNT
          Z[j][i] += Z[k - 1][i] * ZB[j][k];
        #endif          
      }
    }
  }
}

void solveZB(int i, int j, dcomplex x, char sequence[MAXSIZE], int *basePairs, int **basePairCounts, dcomplex **ZB, dcomplex **ZM) { 
  // (i, j) assumed to b.p. in here.
  int k, l;
  
  // In a hairpin, (i + 1, j - 1) all unpaired.
  ZB[i][j] += pow(x, basePairCounts[i][j] + jPairedTo(i, j, basePairs));
  // ZB[i][j] += exp(-HP_Energy(i, j, S0, sequence + 1) / kT) * pow(x, basePairCounts[i][j] + jPairedTo(i, j, basePairs));
  #ifdef STRUCTURE_COUNT
    ZB[j][i] += 1;
  #endif
  
  // Interior loop / bulge / stack / multiloop.
  for (k = i + 1; k <= j - MIN_PAIR_DIST - 1; ++k) {
    for (l = MAX(k + MIN_PAIR_DIST + 1, j - MAX_INTERIOR_DIST - 1); l < j; ++l) {
      if (BP(k, l, sequence)) {
        // In interior loop / bulge / stack with (i, j) and (k, l), (i + 1, k - 1) and (l + 1, j - 1)
        // are all unpaired.
        ZB[i][j] += ZB[k][l] * pow(x, basePairCounts[i][j] - basePairCounts[k][l] + jPairedTo(i, j, basePairs));
        // ZB[i][j] += ZB[k][l] * exp(-IL_Energy(i, j, k, l, S0) / kT) * pow(x, basePairCounts[i][j] - basePairCounts[k][l] + jPairedTo(i, j, basePairs));
        #ifdef STRUCTURE_COUNT
          ZB[j][i] += ZB[l][k];
        #endif
        
        if (k > i + MIN_PAIR_DIST + 2) {
          // If (i, j) is the closing b.p. of a multiloop, and (k, l) is the rightmost base pair, 
          // there is at least one hairpin between (i + 1, k - 1).
          ZB[i][j] += ZM[i + 1][k - 1] * ZB[k][l] * pow(x, basePairCounts[i][j] - basePairCounts[i + 1][k - 1] - basePairCounts[k][l] + jPairedTo(i, j, basePairs));
          // ZB[i][j] += exp(-(ML_close + MLbasepairAndAUpenalty(j, i, S0)) / kT) * ZB[k][l] * ZM[i + 1][k - 1] * pow(x, basePairCounts[i][j] - basePairCounts[i + 1][k - 1] - basePairCounts[k][l] + jPairedTo(i, j, basePairs));
          #ifdef STRUCTURE_COUNT
            ZB[j][i] += ZM[k - 1][i + 1] * ZB[l][k];
          #endif
        }
      }
    }
  }
}

void solveZM(int i, int j, dcomplex x, char sequence[MAXSIZE], int *basePairs, int **basePairCounts, dcomplex **ZB, dcomplex **ZM) { 
  int k;
  
  ZM[i][j] += ZM[i][j - 1] * pow(x, jPairedIn(i, j, basePairs));
  // ZM[i][j] += ZM[i][j - 1] * exp(-1 / kT) * pow(x, jPairedIn(i, j, basePairs));
  
  #ifdef STRUCTURE_COUNT
    ZM[j][i] += ZM[j - 1][i];
  #endif
  
  for (k = i; k <= j - MIN_PAIR_DIST - 1; ++k) {
    if (BP(k, j, sequence)) {
      // Only one stem.
      ZM[i][j] += ZB[k][j] * pow(x, basePairCounts[i][j] - basePairCounts[k][j]);
      // ZM[i][j] += ZB[k][j] * exp(-ML_base * (k - i) / kT) * pow(x, basePairCounts[i][j] - basePairCounts[k][j]);
      
      #ifdef STRUCTURE_COUNT
        ZM[j][i] += ZB[j][k];
      #endif
      
      // k needs to be greater than MIN_PAIR_DIST + 2 from i to fit more than one stem.
      if (k > i + MIN_PAIR_DIST + 2) {
        ZM[i][j] += ZM[i][k - 1] * ZB[k][j] * pow(x, basePairCounts[i][j] - basePairCounts[i][k - 1] - basePairCounts[k][j]);
        // ZM[i][j] += ZM[i][k - 1] * ZB[k][j] * exp(-ML_base / kT) * pow(x, basePairCounts[i][j] - basePairCounts[i][k - 1] - basePairCounts[k][j]);
        
        #ifdef STRUCTURE_COUNT
          ZM[j][i] += ZM[k - 1][i] * ZB[j][k];
        #endif
      }
    }
  }
}

void solveLinearSystem(dcomplex **rootsOfUnity) {
  int i, j;
  
  std::cout << "Roots:" << std::endl;
  
  for (i = 0; i <= seqlen; ++i) {
    std::cout << i << ": (" << rootsOfUnity[i][0].real() << ", i * " << rootsOfUnity[i][0].imag() << ") -> (" << rootsOfUnity[i][1].real() << ", i * " << rootsOfUnity[i][1].imag() << ")" << std::endl;
  }
  
  std::cout << std::endl << std::endl;
  
  // Might need to free this memory.
  LaGenMatComplex A(seqlen + 1, seqlen + 1);
  LaVectorComplex X(seqlen + 1);
  LaVectorComplex B(seqlen + 1);
  
  for (i = 0; i <= seqlen; ++i) {
    for (j = 0; j <= seqlen; ++j) {
      A(i, j).r = pow(rootsOfUnity[i][0], j).real();
      A(i, j).i = pow(rootsOfUnity[i][0], j).imag();
    }
    
    B(i).r = rootsOfUnity[i][1].real();
    B(i).i = rootsOfUnity[i][1].imag();
  }

  // std::cout << "Before:" << std::endl;
  // std::cout << A << std::endl << std::endl << std::endl;
  // std::cout << X << std::endl << std::endl << std::endl;
  // std::cout << B << std::endl << std::endl << std::endl;
  
  LaLinearSolveIP(A, X, B);
  
  std::cout << "Solution:" << std::endl;
  
  for (i = 0; i <= seqlen; ++i) {
    std::cout << i << ": " << X(i) << std::endl;
  }
  // std::cout << A << std::endl << std::endl << std::endl;
  // std::cout << X << std::endl << std::endl << std::endl;
  // std::cout << B << std::endl << std::endl << std::endl;
}

int jPairedTo(int i, int j, int *basePairs) {
  return basePairs[i] == j ? -1 : 1;
}

int jPairedIn(int i, int j, int *basePairs) {
  return basePairs[j] >= i && basePairs[j] < j ? 1 : 0;
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
