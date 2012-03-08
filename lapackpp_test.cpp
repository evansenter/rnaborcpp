#define LA_COMPLEX_SUPPORT 1

#include <stdio.h>
#include <lapackpp.h>
#include <iostream>

using namespace std;

typedef complex<double> dcomplex;

int main(int argc, char* argv[]) {
  int row = 1000;
  int col = 1000;
  
  dcomplex **squareMatrix = new dcomplex*[row];
  for (int i = 0; i < col; i++) {
    squareMatrix[i] = new dcomplex[col];
  }  
    
  LaGenMatComplex A(row, col);
  LaVectorComplex X(col);
  LaVectorComplex B(col);
  
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      squareMatrix[i][j] = dcomplex(i, j);
      
      A(i, j).r = i + j + 0.5;
      A(i, j).i = -(i + j + 0.5);
    }
    
    B(i).r = i + 0.5;
    B(i).i = -(i + 0.5);
  }
  // 
  // std::cout << A << std::endl;
  // 
  // std::complex<double> *x = new std::complex<double>(3, 2);
  // std::complex<double> *y = new std::complex<double>(4, 5);
  // 
  // x + y;
  // 
  // std::cout << x -> real() << std::endl;
  
  LaLinearSolveIP(A, X, B);
  
  cout << X << endl;
	  
  return 0;
}

// c++ lapackpp_test.cpp -I /usr/local/include/lapackpp `pkg-config lapackpp --libs` -o lapackpp_test && ./lapackpp_test