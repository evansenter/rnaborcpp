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
#define SCALING_FACTOR 2.5
#define MIN_PAIR_DIST 3
#define MAX_INTERIOR_DIST 30
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define ZERO_C dcomplex(0.0, 0.0)
#define ONE_C dcomplex(1.0, 0.0)
#define DEBUG 1

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
        Z[i][i + d] = dcomplex(1 / pow(SCALING_FACTOR, d + 1), 0);
          
        if (STRUCTURE_COUNT && d > 0) {
          Z[i + d][i] = ONE_C;
        }
      }
    }
    
    if (DEBUG && root == 0) {
      printMatrix(Z, (char *)"Initialized matrix (1-indexed):", 1, seqlen, 1, seqlen);
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
    
    if (DEBUG && root == 0) {
      printMatrix(Z, (char *)"Evaluated matrix (1-indexed, zeroth root):", 1, seqlen, 1, seqlen);
    }
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
  
  Z[i][j] += Z[i][j - 1] * pow(x, jPairedIn(i, j, basePairs)) / dcomplex(pow(SCALING_FACTOR, 1), 0);
    
  if (STRUCTURE_COUNT) {
    // Conditional if j - 1 == i
    Z[j][i] += (j - 1 == i ? ONE_C : Z[j - 1][i]);
  }
    
  for (k = i; k <= j - MIN_PAIR_DIST - 1; ++k) { 
    // (k, j) is the rightmost base pair in (i, j).
    if (BP(k, j, sequence)) {
      if (k == i) {
        Z[i][j] += ZB[k][j] * exp(-AU_Penalty(i, j, S0) / kT) * pow(x, jPairedIn(i, j, basePairs));
          
        if (STRUCTURE_COUNT) {
          Z[j][i] += ZB[j][k];
        }
      } else {
        Z[i][j] += Z[i][k - 1] * ZB[k][j] * exp(-AU_Penalty(k, j, S0) / kT) * pow(x, basePairCounts[i][j] - basePairCounts[i][k - 1] - basePairCounts[k][j]);
          
        if (STRUCTURE_COUNT) {
          Z[j][i] += Z[k - 1][i] * ZB[j][k];
        }
      }
    }
  }
}

void solveZB(int i, int j, dcomplex x, char sequence[MAXSIZE], int *basePairs, int **basePairCounts, dcomplex **ZB, dcomplex **ZM) { 
  // (i, j) assumed to b.p. in here.
  int k, l;
  
  // In a hairpin, (i + 1, j - 1) all unpaired.
  ZB[i][j] += exp(-HP_Energy(i, j, S0, sequence + 1) / kT) * pow(x, basePairCounts[i][j] + jPairedTo(i, j, basePairs));
  
  if (STRUCTURE_COUNT) {
    ZB[j][i] += 1;
  }
  
  // Interior loop / bulge / stack / multiloop.
  for (k = i + 1; k <= j - MIN_PAIR_DIST - 1; ++k) {
    for (l = MAX(k + MIN_PAIR_DIST + 1, j - MAX_INTERIOR_DIST - 1); l < j; ++l) {
      if (BP(k, l, sequence)) {
        // In interior loop / bulge / stack with (i, j) and (k, l), (i + 1, k - 1) and (l + 1, j - 1)
        // are all unpaired.
        ZB[i][j] += ZB[k][l] * exp(-IL_Energy(i, j, k, l, S0) / kT) * pow(x, basePairCounts[i][j] - basePairCounts[k][l] + jPairedTo(i, j, basePairs)) / dcomplex(pow(SCALING_FACTOR, j - l + k - i), 0);
        
        if (STRUCTURE_COUNT) {
          ZB[j][i] += ZB[l][k];
        }
        
        if (k > i + MIN_PAIR_DIST + 2) {
          // If (i, j) is the closing b.p. of a multiloop, and (k, l) is the rightmost base pair, 
          // there is at least one hairpin between (i + 1, k - 1).
          ZB[i][j] += exp(-(ML_close + MLbasepairAndAUpenalty(j, i, S0)) / kT) * ZM[i + 1][k - 1] * ZB[k][l] * pow(x, basePairCounts[i][j] - basePairCounts[i + 1][k - 1] - basePairCounts[k][l] + jPairedTo(i, j, basePairs)) / dcomplex(pow(SCALING_FACTOR, j - l + 1), 0);
          
          if (STRUCTURE_COUNT) {
            ZB[j][i] += ZM[k - 1][i + 1] * ZB[l][k];
          }
        }
      }
    }
  }
}

void solveZM(int i, int j, dcomplex x, char sequence[MAXSIZE], int *basePairs, int **basePairCounts, dcomplex **ZB, dcomplex **ZM) { 
  int k;
  
  ZM[i][j] += ZM[i][j - 1] * pow(x, jPairedIn(i, j, basePairs)) / dcomplex(pow(SCALING_FACTOR, 1), 0);
  // ZM[i][j] += ZM[i][j - 1] * exp(-1 / kT) * pow(x, jPairedIn(i, j, basePairs)) / dcomplex(pow(SCALING_FACTOR, 1), 0);
  
  if (STRUCTURE_COUNT) {
    ZM[j][i] += ZM[j - 1][i];
  }
  
  for (k = i; k <= j - MIN_PAIR_DIST - 1; ++k) {
    if (BP(k, j, sequence)) {
      // Only one stem.
      ZM[i][j] += ZB[k][j] * exp(-ML_base * (k - i) / kT) * pow(x, basePairCounts[i][j] - basePairCounts[k][j]) / dcomplex(pow(SCALING_FACTOR, k - i), 0);
      
      if (STRUCTURE_COUNT) {
        ZM[j][i] += ZB[j][k];
      }
      
      // k needs to be greater than MIN_PAIR_DIST + 2 from i to fit more than one stem.
      if (k > i + MIN_PAIR_DIST + 2) {
        ZM[i][j] += ZM[i][k - 1] * ZB[k][j] * exp(-ML_base / kT) * pow(x, basePairCounts[i][j] - basePairCounts[i][k - 1] - basePairCounts[k][j]);
        
        if (STRUCTURE_COUNT) {
          ZM[j][i] += ZM[k - 1][i] * ZB[j][k];
        }
      }
    }
  }
}

void solveLinearSystem(dcomplex **rootsOfUnity) {
  int i, j;
  dcomplex sum = ZERO_C, poweredRoot;
  
  if (DEBUG) {
    printMatrix(rootsOfUnity, (char *)"Roots and solutions:", 0, seqlen, 0, 1);
  }
  
  // Might need to free this memory.
  LaGenMatComplex A(seqlen + 1, seqlen + 1);
  LaVectorComplex X(seqlen + 1);
  LaVectorComplex B(seqlen + 1);
  
  for (i = 0; i <= seqlen; ++i) {
    for (j = 0; j <= seqlen; ++j) {
      poweredRoot = pow(rootsOfUnity[i][0], j);
      
      A(i, j).r = poweredRoot.real();
      A(i, j).i = poweredRoot.imag();
    }
    
    B(i).r = rootsOfUnity[i][1].real();
    B(i).i = rootsOfUnity[i][1].imag();
  }

  // std::cout << "Before:" << std::endl;
  // std::cout << A << std::endl << std::endl << std::endl;
  // std::cout << X << std::endl << std::endl << std::endl;
  // std::cout << B << std::endl << std::endl << std::endl;
  
  LaLinearSolveIP(A, X, B);

  for (i = 0; i <= seqlen; ++i) {
    sum = sum + dcomplex(X(i).r, X(i).i);
  }
  
  std::cout << "Solution:" << std::endl;
  std::cout << "Sum: " << sum << std::endl;
  
  for (i = 0; i <= seqlen; ++i) {
    std::cout << i << ": " << dcomplex(X(i).r, X(i).i) / sum << std::endl;
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

void printMatrix(dcomplex **matrix, char *title, int iStart, int iStop, int jStart, int jStop) {
  int i, j;
  
  printf("%s\n", title);
      
  for (i = iStart; i <= iStop; ++i) {
    for (j = jStart; j <= jStop; ++j) {
      printf("%+.3f, %-+10.3f", matrix[i][j].real(), matrix[i][j].imag());
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}