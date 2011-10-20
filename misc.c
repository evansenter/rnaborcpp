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
#include "misc.h"
#include <limits.h>

/*The following function checks if the ith and jth positions in the sequence could base pair*/
int BP(int i,int j,char sequence[MAXSIZE])  //safe
{if(j-i<4)
    return 0;
  else if ((sequence[i]=='A'&& sequence[j]=='U')||(sequence[i]=='U'&& sequence[j]=='A'))
    return 1;
  else if ((sequence[i]=='C'&& sequence[j]=='G')||(sequence[i]=='G'&& sequence[j]=='C'))
    return 1;
  else if ((sequence[i]=='U'&& sequence[j]=='G')||(sequence[i]=='G'&& sequence[j]=='U'))
    return 1;
  else
    return 0;
}

int ElementInList(char subshape[SHAPELENGTH],char shapeList[SHAPELENGTH][SHAPELENGTH],int listlen) //check if a subshape is in list, if yes,return the index of the subshape,if not return -1
{int i;
 for(i=0;i<listlen;++i)
    {if(strcmp(shapeList[i],subshape)==0)
        return i;
    }
 return -1;
}

int NumInList(int a, int L[SHAPELENGTH],int listlen) //check if a number is in this list,if not return -1
{int i;
 for(i=0;i<listlen;++i)
    {if(a==L[i])
        return i;
    }
 return -1;
}

/*This function returns the largest value in an array.*/
double MaxInArray(double array[],int arraylen)   //safe
{ double largest=array[0];
  int i;
  for (i=0;i<arraylen;++i)
    {if (array[i]>largest)
        largest=array[i];
    }
  return largest;
}

/*The following function check if the input string is a valid RNA sequence*/
int CheckSequence(char sequence[MAXSIZE]){     //safe
  int i;
  for (i=1;i<strlen(sequence);++i)
    {
     if (toupper(sequence[i])!='A' && toupper(sequence[i])!='U' && toupper(sequence[i])!='C' && toupper(sequence[i])!='G'&& toupper(sequence[i])!='T') //check if there are invalid characters
        {
         printf("The input string should only contain A,U,T,C,G!\n");
         exit(1);
        }
     else if (toupper(sequence[i])=='T') //change T to U
        sequence[i]='U';
     else
        sequence[i]=toupper(sequence[i]); // change lower case to upper case
    }
  return 0;
}

double ***Allocate3DMatrix(int a, int b, int c)
// allocate a matrix of size a*b*c
{ int i;
  int j;
  double ***Matrix;
  Matrix=(double ***) malloc((a) * sizeof (double **));
  if(Matrix == NULL)
        {printf("out of memory\n");
         exit(1);}
  for (i=0; i<a; i++){
    Matrix[i]=(double **) malloc((b)*sizeof(double *));
    if(Matrix[i] == NULL)
        {printf("out of memory\n");
         exit(1);}
  }
  for(i=0;i<a;i++)
    {for(j=0;j<b;j++)
        {
        Matrix[i][j]=(double *) malloc((c)*sizeof(double));
        if(Matrix[i][j] == NULL)
        {printf("out of memory\n");
         exit(1);}
        }
    }
  return Matrix;
}

int Free3DMatrix(double ***Matrix, int a, int b, int c)
// Free a matrix of size a*b*c
{ int i;
  int j;
  for(i=0;i<a;i++)
	{
	for(j=0;j<b;j++)
	free(Matrix[i][j]);
	}
  for(i=0;i<a;i++)
	free(Matrix[i]);
  free(Matrix);
  return 0;
}

double **Allocate2DMatrix(int a, int b)
// allocate a matrix of size a*b
{ int i;
  double **Matrix;
  Matrix=(double **) malloc((a) * sizeof (double *));
  if(Matrix == NULL)
        {printf("out of memory\n");
         exit(1);}
  for (i=0; i<a; i++){
    Matrix[i]=(double *) malloc((b)*sizeof(double));
    if(Matrix[i] == NULL)
        {printf("out of memory\n");
         exit(1);}
  }
  return Matrix;
}

int Free2DMatrix(double **Matrix,int a, int b)
// free a matrix of size a*b
{ int i;
  for (i=0; i<a; i++){
    free(Matrix[i]);
  }
  free(Matrix);
  return 0;
}

int pairedIn(int i, int j, int k, char *structure) {
  if (j - i < 4) {
    return 0;
  } else {
    int counter, stackPointer;
    int pairs[j + 1];
    int stack[j - i + 1];

    stackPointer = 0;
    
    for (counter = i; counter <= j; ++counter) {
      pairs[counter] = 0;
    }

    for (counter = i; counter <= j; ++counter) {
      if (structure[counter] == '(') {
        stack[stackPointer++] = counter;
      } else if (structure[counter] == ')' && stackPointer > 0) {
        stackPointer -= 1;
        pairs[stack[stackPointer]] = counter;
        pairs[counter]             = stackPointer;
      }
    }
        
    return pairs[k];
  }
}

void populatePairings(int i, int j, char *structure, int *pairs) {
  int structurePointer, stackPointer;
  int stack[j - i + 1];
  
  stackPointer = 0;
  
  for (structurePointer = i; structurePointer <= j; ++structurePointer) {
    if (structure[structurePointer] == '(') {
      stack[stackPointer++] = structurePointer;
    } else if (structure[structurePointer] == ')' && stackPointer > 0) {
      pairs[stack[--stackPointer]] = structurePointer;
    }
  }
}

int comparePairings(int i, int j, int *pairs, int *comparisonPairs) {
  int counter, pairings;
  
  pairings = 0;
  
  for (counter = i; counter <= j; ++counter) {
    if (pairs[counter] && pairs[counter] != comparisonPairs[counter]) {
      pairings += 1;
    }
  }
  
  for (counter = i; counter <= j; ++counter) {
    if (comparisonPairs[counter] && comparisonPairs[counter] != pairs[counter]) {
      pairings += 1;
    }
  }
  
  return pairings;
}

int bpDifference(int i, int j, int k, int l, char *structure, int bounds[][2], int numBounds) {
  int counter;
  int pairs[j + 1];
  int comparisonPairs[j + 1];
  
  for (counter = i; counter <= j; ++counter) {
    pairs[counter]           = 0;
    comparisonPairs[counter] = 0;
  }
  
  populatePairings(i, j, structure, pairs);
  
  for (counter = 0; counter < numBounds; ++counter) {
    populatePairings(bounds[counter][0], bounds[counter][1], structure, comparisonPairs);
  }
  
  if (k && l) {
    comparisonPairs[k] = l;
  }
  
  return comparePairings(i, j, pairs, comparisonPairs);
}
