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
#define SET_Z(i, j, value) \
Z[i][j] = value; \
Z[j][i] = value;

typedef std::complex<double> dcomplex;

double** runMcCaskill(char sequence[MAXSIZE]) {
  int root, i, j, d;
  double **Z;
  double **ZB;
  double **ZM;
  
  Z  = Allocate2DMatrix(seqlen + 1, seqlen + 1);
  ZB = Allocate2DMatrix(seqlen + 1, seqlen + 1);
  ZM = Allocate2DMatrix(seqlen + 1, seqlen + 1);
  
  dcomplex **cZ           = new dcomplex*[seqlen + 1];
  dcomplex **cZB          = new dcomplex*[seqlen + 1];
  dcomplex **cZM          = new dcomplex*[seqlen + 1];
  dcomplex **rootsOfUnity = new dcomplex*[seqlen + 1];
  
  for (i = 0; i <= seqlen; i++) {
    cZ[i]              = new dcomplex[seqlen + 1];
    cZB[i]             = new dcomplex[seqlen + 1];
    cZM[i]             = new dcomplex[seqlen + 1];
    rootsOfUnity[i]    = new dcomplex[2];
    rootsOfUnity[i][0] = dcomplex(cos(2 * M_PI * i / (seqlen + 1)), sin(2 * M_PI * i / (seqlen + 1)));
  }
  
  for (root = 0; root <= seqlen; ++root) {
    // Flush the matrices.
    for (i = 0; i <= seqlen; ++i) {
      for (j = 0; j <= seqlen; ++j) {
        Z[i][j]  = 0;
        ZB[i][j] = 0;
        ZM[i][j] = 0;
      }  
    }
  
    // Set base case for Z.
    for (d = 0; d <= MIN_PAIR_DIST; ++d) {
      for (i = 1; i <= seqlen - d; ++i) {
        SET_Z(i, i + d, 1)
      }
    }
    
    for (d = MIN_PAIR_DIST + 1; d < seqlen; ++d) {
      for (i = 1; i <= seqlen - d; ++i) {
        j = i + d;
      
        if (BP(i, j, sequence)) {
          solveZB(i, j, sequence, ZB, ZM);
        }
        
        solveZM(i, j, sequence, ZB, ZM);
      
        solveZ(i, j, sequence, Z, ZB);
      }
    }
  
    // LaGenMatComplex A(row, col);
    // LaVectorComplex X(col);
    // LaVectorComplex B(col);
    // 
    // for (i = 0; i < row; i++) {
    //   for (j = 0; j < col; j++) {
    //     squareMatrix[i][j] = dcomplex(i, j);
    //     
    //     A(i, j).r = i + j + 0.5;
    //     A(i, j).i = -(i + j + 0.5);
    //   }
    //   
    //   B(i).r = i + 0.5;
    //   B(i).i = -(i + 0.5);
    // }
    // 
    // std::cout << A << std::endl;
    // 
    // std::complex<double> *x = new std::complex<double>(3, 2);
    // std::complex<double> *y = new std::complex<double>(4, 5);
    // 
    // x + y;
    // 
    // std::cout << x -> real() << std::endl;
    // 
    // LaLinearSolveIP(A, X, B);
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

void solveZB(int i, int j, char sequence[MAXSIZE], double **ZB, double **ZM) { 
  // (i, j) assumed to b.p. in here.
  int k, l;
  
  // In a hairpin, (i + 1, j - 1) all unpaired.
  ZB[i][j] += exp(-HP_Energy(i, j, S0, sequence + 1) / kT);
  ZB[j][i] += 1;
  
  // Interior loop / bulge / stack / multiloop.
  for (k = i + 1; k <= j - MIN_PAIR_DIST - 1; ++k) {
    for (l = MAX(k + MIN_PAIR_DIST + 1, j - MAX_INTERIOR_DIST - 1); l < j; ++l) {
      if (BP(k, l, sequence)) {
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
    if (BP(k, j, sequence)) {
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