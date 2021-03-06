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
#include "RNAbor.h"
#include "misc.h"
#include "McCaskill.h"
#include <fftw3.h>
#include <iostream>
#define STRUCTURE_COUNT 1
#define MIN_PAIR_DIST 3
#define MAX_INTERIOR_DIST 30
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define ZERO_C dcomplex(0.0, 0.0)
#define ONE_C dcomplex(1.0, 0.0)
#define SCALE(power) pow(SCALING_FACTOR, power)
#define DEBUG 1
#define PRINT_MATRICES 0
#define PRINT_DETAILED_MATRICIES 0
#define FFTW_REAL 0
#define FFTW_IMAG 1

dcomplex** runMcCaskill(char sequence[MAXSIZE], char structure[MAXSIZE]) {
  // Variable declarations.
  int root, i, j, k, d, *basePairs, **bpCounts;
  double scalingFactor;
  
  dcomplex **Z            = new dcomplex*[seqlen + 1];
  dcomplex **ZB           = new dcomplex*[seqlen + 1];
  dcomplex **ZM           = new dcomplex*[seqlen + 1];
  dcomplex **rootsOfUnity = new dcomplex*[seqlen + 1];
  double *coefficients    = new double[seqlen + 1];
  
  // Matrix allocation.
  for (i = 0; i <= seqlen; ++i) {
    Z[i]               = new dcomplex[seqlen + 1];
    ZB[i]              = new dcomplex[seqlen + 1];
    ZM[i]              = new dcomplex[seqlen + 1];
    rootsOfUnity[i]    = new dcomplex[2];
    rootsOfUnity[i][0] = dcomplex(cos(2 * M_PI * i / (seqlen + 1)), sin(2 * M_PI * i / (seqlen + 1)));
  }
  
  if (DEBUG) {
    std::cout << "Structure: ";
    for (i = 0; i < seqlen; ++i) {
      std::cout << structure[i];
    }
    std::cout << std::endl;
  }
  
  basePairs = getBasePairList(structure);
  bpCounts  = fillBasePairCounts(basePairs, seqlen);
  
  if (PRINT_DETAILED_MATRICIES) {
    std::cout << "Base pairs array:" << std::endl;
    std::cout << "[ ";
    for (i = 1; i <= seqlen; ++i) {
      std::cout << basePairs[i] << " ";
    }
    std::cout << "]" << std::endl;
  
    std::cout << "Base pair counts:" << std::endl;
    for (i = 1; i <= seqlen; ++i) {
      std::cout << "[ ";
      for (j = 1; j <= seqlen; ++j) {
        std::cout << (j < i ? 0 : bpCounts[i][j]) << " ";
      }
      std::cout << "]" << std::endl;
    }
  }
  
  // Start main recursions (root <= round(seqlen / 2.0) is an optimization for roots of unity).
  for (root = 0; root <= round(seqlen / 2.0); ++root) {
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
        j = i + d;
        
        Z[i][j] = SCALE(-d);
        
        if (STRUCTURE_COUNT && i != j) {
          Z[j][i] = ONE_C;
        }
      }
    }
    
    if (PRINT_MATRICES) {
      std::cout << "SOLUTIONS FOR ROOT " << root << ": " << rootsOfUnity[root][0] << std::endl;
      
      printMatrix(Z, (char *)"Initialized matrix (1-indexed):", 0, seqlen, 0, seqlen);
    }
    
    for (d = MIN_PAIR_DIST + 1; d < seqlen; ++d) {
      for (i = 1; i <= seqlen - d; ++i) {
        j = i + d;
      
        if (BP(i, j, sequence)) {
          solveZB(i, j, rootsOfUnity[root][0], sequence, basePairs, bpCounts, ZB, ZM);
        }
        
        solveZM(i, j, rootsOfUnity[root][0], sequence, basePairs, bpCounts, ZB, ZM);
        
        solveZ(i, j, rootsOfUnity[root][0], sequence, basePairs, bpCounts, Z, ZB);
      }
    }
    
    rootsOfUnity[root][1] = Z[1][seqlen];
    
    if (!root) {
      scalingFactor = Z[1][seqlen].real() * SCALE(seqlen - 1);
    }
    
    if (PRINT_MATRICES) {
      printMatrix(Z, (char *)"Evaluated matrix (1-indexed, zeroth root):", 0, seqlen, 0, seqlen);
    }
    
    if (DEBUG && root == 0) {
      printf("c:         %f\n", SCALING_FACTOR);
      printf("c^(n-1):   %f\n", SCALE(seqlen - 1));
      printf("n:         %d\n", seqlen);
      printf("x:         %+f, %+f\n", rootsOfUnity[root][0].real(), rootsOfUnity[root][0].imag());
      printf("Q[1][%d]:  %.15f\n", seqlen, Z[1][seqlen].real());
      printf("Z[1][%d]:  %.15f\n", seqlen, Z[1][seqlen].real() * SCALE(seqlen - 1));
      printf("Z[%d][1]:  %.15f\n\n", seqlen, Z[seqlen][1].real());
      printf("Working");
    }
    
    if (DEBUG) {
      std::cout << '.' << std::flush;
    }
  }
  
  // Optimization leveraging complementarity of roots of unity.
  if (seqlen % 2) {
    i = root - 2;
  } else {
    i = root - 1;
  }
  
  for (; root <= seqlen && i > 0; --i, ++root) {
    rootsOfUnity[root][1] = dcomplex(rootsOfUnity[i][1].real(), -rootsOfUnity[i][1].imag());
  }
  
  if (DEBUG) {
    printf("\n\n");
  }
  
  solveSystem(rootsOfUnity, coefficients, scalingFactor);
  
  return Z;
}

void solveZ(int i, int j, dcomplex x, char sequence[MAXSIZE], int *basePairs, int **bpCounts, dcomplex **Z, dcomplex **ZB) { 
  int k, delta;
  
  delta    = jPairedIn(i, j, basePairs);
  Z[i][j] += (Z[i][j - 1] * pow(x, delta)) / SCALE(1);
    
  if (STRUCTURE_COUNT) {
    Z[j][i] += Z[j - 1][i];
  }
    
  for (k = i; k <= j - MIN_PAIR_DIST - 1; ++k) { 
    // (k, j) is the rightmost base pair in (i, j).
    if (BP(k, j, sequence)) {
      if (k == i) {
        delta    = bpCounts[i][j] - bpCounts[k][j];
        Z[i][j] += ZB[k][j] * pow(x, delta) * exp(-AU_Penalty(i, j, S0) / kT);
          
        if (STRUCTURE_COUNT) {
          Z[j][i] += ZB[j][k];
        }
      } else {
        delta    = bpCounts[i][j] - bpCounts[i][k - 1] - bpCounts[k][j];
        Z[i][j] += (Z[i][k - 1] * ZB[k][j] * pow(x, delta) * exp(-AU_Penalty(k, j, S0) / kT)) / SCALE(1);
          
        if (STRUCTURE_COUNT) {
          Z[j][i] += Z[k - 1][i] * ZB[j][k];
        }
      }
    }
  }
}

void solveZB(int i, int j, dcomplex x, char sequence[MAXSIZE], int *basePairs, int **bpCounts, dcomplex **ZB, dcomplex **ZM) { 
  // (i, j) assumed to b.p. in here.
  int k, l, delta;
  
  // In a hairpin, [i + 1, j - 1] unpaired.
  delta     = bpCounts[i][j] + jPairedTo(i, j, basePairs);
  ZB[i][j] += (pow(x, delta) * exp(-HP_Energy(i, j, S0, sequence + 1) / kT)) / SCALE(j - i);
  
  if (STRUCTURE_COUNT) {
    ZB[j][i] += 1;
  }
  
  // Interior loop / bulge / stack / multiloop.
  for (k = i + 1; k <= j - MIN_PAIR_DIST - 2; ++k) {
    for (l = k + MIN_PAIR_DIST + 1; l < j; ++l) {
      if (BP(k, l, sequence)) {
        if (k - i + j - l - 2 <= MAX_INTERIOR_DIST) {
          // In interior loop / bulge / stack with (i, j) and (k, l), (i + 1, k - 1) and (l + 1, j - 1)
          // are all unpaired.
          delta     = bpCounts[i][j] - bpCounts[k][l] + jPairedTo(i, j, basePairs);
          ZB[i][j] += (ZB[k][l] * pow(x, delta) * exp(-IL_Energy(i, j, k, l, S0) / kT)) / SCALE(j - l + k - i);
        
          if (STRUCTURE_COUNT) {
            ZB[j][i] += ZB[l][k];
          }
        }
        
        if (k > i + MIN_PAIR_DIST + 2) {
          // If (i, j) is the closing b.p. of a multiloop, and (k, l) is the rightmost base pair, 
          // there is at least one hairpin between (i + 1, k - 1).
          delta     = bpCounts[i][j] - bpCounts[i + 1][k - 1] - bpCounts[k][l] + jPairedTo(i, j, basePairs);
          ZB[i][j] += (ZM[i + 1][k - 1] * ZB[k][l] * pow(x, delta) * exp(-(ML_close + MLbasepairAndAUpenalty(j, i, S0)) / kT)) / SCALE(j - l + 2);
          
          if (STRUCTURE_COUNT) {
            ZB[j][i] += ZM[k - 1][i + 1] * ZB[l][k];
          }
        }
      }
    }
  }
}

void solveZM(int i, int j, dcomplex x, char sequence[MAXSIZE], int *basePairs, int **bpCounts, dcomplex **ZB, dcomplex **ZM) { 
  int k, delta;
  
  delta     = jPairedIn(i, j, basePairs);
  ZM[i][j] += (ZM[i][j - 1] * pow(x, delta) * exp(-1 / kT)) / SCALE(1);
  
  if (STRUCTURE_COUNT) {
    ZM[j][i] += ZM[j - 1][i];
  }
  
  for (k = i; k <= j - MIN_PAIR_DIST - 1; ++k) {
    if (BP(k, j, sequence)) {
      // Only one stem.
      delta     = bpCounts[i][j] - bpCounts[k][j];
      ZM[i][j] += (ZB[k][j] * pow(x, delta) * exp(-ML_base * (k - i) / kT)) / SCALE(k - i);
      
      if (STRUCTURE_COUNT) {
        ZM[j][i] += ZB[j][k];
      }
      
      // More than one stem.
      if (k > i) {
        delta     = bpCounts[i][j] - bpCounts[i][k - 1] - bpCounts[k][j];
        ZM[i][j] += (ZM[i][k - 1] * ZB[k][j] * pow(x, delta) * exp(-ML_base / kT)) / SCALE(1);
        
        if (STRUCTURE_COUNT) {
          ZM[j][i] += ZM[k - 1][i] * ZB[j][k];
        }
      }
    }
  }
}

void solveSystem(dcomplex **rootsOfUnity, double *coefficients, double scalingFactor) {
  int i;
  dcomplex sum = ZERO_C;
  
  if (DEBUG) {
    printMatrix(rootsOfUnity, (char *)"START ROOTS AND SOLUTIONS", 0, seqlen, 0, 1);
    std::cout << "END ROOTS AND SOLUTIONS" << std::endl << std::endl;
    std::cout << "Scaling factor (Z{1, n}): " << scalingFactor << std::endl;
  }

  fftw_complex signal[seqlen + 1];
  fftw_complex result[seqlen + 1];
  
  for (i = 0; i <= seqlen; i++) {
    signal[i][FFTW_REAL] = (pow(10, PRECISION) * rootsOfUnity[i][1].real()) / scalingFactor;
    signal[i][FFTW_IMAG] = (pow(10, PRECISION) * rootsOfUnity[i][1].imag()) / scalingFactor;
  }
  
  fftw_plan plan = fftw_plan_dft_1d(seqlen + 1, signal, result, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  
  std::cout << "START DISTRIBUTION" << std::endl;
  for (i = 0; i <= seqlen; i++) {
    coefficients[i] = PRECISION == 0 ? result[i][FFTW_REAL] / (seqlen + 1) : pow(10.0, -PRECISION) * static_cast<int>(result[i][FFTW_REAL] / (seqlen + 1));
    sum            += coefficients[i];
    
    std::cout << i << ": " << coefficients[i] << std::endl;
  }
  std::cout << "END DISTRIBUTION" << std::endl << std::endl;
  
  std::cout << "Sum: " << sum << std::endl;
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
  int i, d, **bpCounts;
  
  bpCounts = (int **) calloc(n + 1, sizeof(int *));
  for (i = 1; i <= n; ++i) {
    bpCounts[i] = (int *) calloc(n + 1, sizeof(int));
  }
  
  for (d = MIN_PAIR_DIST + 1; d < n; ++d) {
    for (i = 1; i <= n - d; ++i) {
      bpCounts[i][i + d] = numberOfBasePairs(i, i + d, basePairs);
    }
  }
  
  return bpCounts;
}

void printMatrix(dcomplex **matrix, char *title, int iStart, int iStop, int jStart, int jStop) {
  int i, j;
  
  printf("%s\n", title);
      
  for (i = iStart; i <= iStop; ++i) {
    for (j = jStart; j <= jStop; ++j) {
      printf("%+.15f, %-+25.15f", matrix[i][j].real(), matrix[i][j].imag());
    }
    std::cout << std::endl;
  }
}
